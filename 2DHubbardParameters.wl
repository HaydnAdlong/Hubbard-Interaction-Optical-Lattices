(* ::Package:: *)

(* ::Title:: *)
(*Lattice2DParameters Package*)


(* ::Subtitle:: *)
(*A package for determining the Hubbard onsite interaction energy U for a quasi-2D optical lattice*)


(* ::Author:: *)
(*Haydn S Adlong*)
(*Package version 0.1*)
(*Jan 2022*)


(* ::Section:: *)
(*Introduction*)


(*BeginPackage[ "Lattice2DParameters`" ]*)

hubbardU::usage = " hubbardU[] determines the onsite interaction U, as a function of  "


(* ::Text:: *)
(*In this package we use units where m = \[HBar] = d = 1 with*)
(**m the particle mass,*)
(**\[HBar] the reduced Planck's constant and*)
(**d the lattice spacing.*)
(**)
(*In these units the recoil energy is *)


Vrec=\[Pi]^2/2.;


(* ::Text:: *)
(*The standard notation includes*)
(**p: incoming relative momentum*)
(**k: relative momentum for integrating*)
(**t\[Sigma]: hopping parameter of spin-\[Sigma] species*)


(* ::Section:: *)
(*1. Useful packages and functions*)


Needs @ "NumericalDifferentialEquationAnalysis`";
Needs @ "NumericalCalculus`";

sortEigensystemByMagnitude[es_] := Sort[es\[Transpose], #1[[1]] < #2[[1]] &]\[Transpose]; (* sorts eigensystems so largest eigenvalue is last *)

TensorMultiply[A_, B_, pairs_] := Activate @ TensorContract[Inactive[TensorProduct][A, B], (# + {0, TensorRank[A]}) &/@ pairs];

intervalInverse[Interval[int___]] := Interval @@ Partition[Flatten @ {int} /. {{-\[Infinity], mid___, \[Infinity]} :> {mid}, {-\[Infinity], mid__} :> {mid, \[Infinity]}, {mid__, \[Infinity]} :> {-\[Infinity], mid}, {mid___} :> {-\[Infinity], mid, \[Infinity]}}, 2]

intervalComplement[a_Interval, b__Interval] := IntervalIntersection[a, intervalInverse @ IntervalUnion[b]]

Options[FindAllCrossings] = Sort[Join[Options[FindRoot], {MaxRecursion -> Automatic, PerformanceGoal :> $PerformanceGoal, PlotPoints -> Automatic}]];

FindAllCrossings[f_, {t_, a_, b_}, opts___] := Module[{r, s, s1, ya},
	{r, ya} = Transpose[First[Cases[Normal[Plot[f, {t, a, b}, Method -> Automatic, Evaluate[Sequence @@ FilterRules[Join[{opts}, Options[FindAllCrossings]], Options[Plot]]]]], Line[l_] :> l, Infinity]]];
	s1 = Sign[ya]; If[!MemberQ[Abs[s1], 1], Return[{}]];
	s = Times @@@ Partition[s1, 2, 1];
	If[MemberQ[s, -1] || MemberQ[Take[s, {2, -2}], 0],
		Union[Join[Pick[r, s1, 0],
		Select[t /. Map[FindRoot[f, {t, r[[#]], r[[#+1]]},
		Evaluate[Sequence @@ FilterRules[Join[{opts}, Options[FindAllCrossings]], Options[FindRoot]]]] &, Flatten[Position[s, -1]]], a <= # <= b &]]], {}
	]
]

quasi2DSum0toInf[var_, \[Omega]z_] := (Sqrt[\[Pi]] Gamma[1+var/(2 \[Omega]z)])/(var Gamma[(var+\[Omega]z)/(2 \[Omega]z)]); (*Sum[(2s)!/(2^s s!)^2 1/(var+2\[Omega]z s),{s,0,\[Infinity]}]*)
quasi2DSum0toInfApprox[var_, \[Omega]z_] := Sqrt[\[Pi]]/var  ((var/(2 \[Omega]z))^2 + var/(4 \[Omega]z) + 1/8)^(1/4);
quasi2DSum1toInf[var_, \[Omega]z_] := (-1+(Sqrt[\[Pi]] Gamma[1+var/(2 \[Omega]z)])/Gamma[(var+\[Omega]z)/(2 \[Omega]z)])/var; (*Sum[(2s)!/(2^s s!)^2 1/(var+2\[Omega]z s),{s,1,\[Infinity]}]*)
quasi2DSum1toInfApprox[var_, \[Omega]z_] := -1/var + Sqrt[\[Pi]]/var  ((var/(2 \[Omega]z))^2 + var/(4 \[Omega]z) + 1/8)^(1/4);


(* ::Section:: *)
(*2. Integration grids*)


(* ::Text:: *)
(*	The grids which will be used to the integral equations.*)
(**)
(*Regular integration from 0 to \[Pi]*)


setupKIntReg[nk_] := Module[{kTemp, dkTemp, kReg, dkReg},
	{kReg, dkReg} = Transpose @ GaussianQuadratureWeights[nk, -\[Pi], \[Pi]]
];


(* ::Text:: *)
(*Principal Value (PV) integration from k0 to k1 evenly distributed around the poles*)


setupKIntPV[nk_, poles_, k0_, k1_] := Module[
	{
		kTemp, dkTemp, poleSpacing, minPoleSpacing, minPoleToEndpointSpacing,
		integrationSpacing, integrationSubDomainsTemp, integrationSubDomains, kPV, dkPV
	},
	poleSpacing = Table[
		Part[poles, i + 1] + -Part[poles, i],
		{i, Length[poles] - 1}
	];
	minPoleSpacing = Min @ Abs @ poleSpacing;
	minPoleToEndpointSpacing = Min[
		Abs[
			{k0 + -Part[poles, 1], k1 + -Part[poles, -1]}
		]
	];
	If[
		Length @ poles == 1,
		integrationSpacing = minPoleToEndpointSpacing;,
		If[
			minPoleSpacing > minPoleToEndpointSpacing,
			integrationSpacing = minPoleSpacing / Ceiling[minPoleSpacing / minPoleToEndpointSpacing];,
			integrationSpacing = minPoleSpacing;
		]
	];

	integrationSubDomainsTemp = DeleteDuplicates[
		Round[
			Sort[
				Flatten[
					Join[
						{
							N @ k0, poles -integrationSpacing, poles + integrationSpacing,
							poles, N @ k1
						}
					]
				]
			],
			10 ^ -9
		]
	];
	integrationSubDomains = N[
		Partition[
			Sort[
				Join[integrationSubDomainsTemp, integrationSubDomainsTemp[[2;;]]]
			],
			2
		]
	];
	{kTemp, dkTemp} = Transpose @ GaussianQuadratureWeights[nk, 0, 1];
	{kPV, dkPV} = {
		Flatten[
			Table[
				Plus[
					Times[kTemp,
						Part[integrationSubDomains, i, 2] + -Part[integrationSubDomains, i, 1]
					],
					integrationSubDomains[[i, 1]]
				],
				{i, Length @ integrationSubDomains}
			]
		],
		Flatten[
			Table[
				Times[dkTemp,
					Part[integrationSubDomains, i, 2] + -Part[integrationSubDomains, i, 1]
				],
				{i, Length @ integrationSubDomains}
			]
		]
	}
];


(* ::Section:: *)
(*3. Hubbard model integration*)


(* ::Text:: *)
(*	Here we perform the integral required to calculated the Hubbard model T matrix. The integral is performed in polar coordinates.*)


HubbardInvFreeGreenFuncPolarCoord[tupx_, tdownx_, tupy_, tdowny_, px_, py_, k_, \[Theta]_] :=
	2 tupx (Cos[-k Cos[\[Theta]]] - Cos[px]) + 2 tupy (Cos[-k Sin[\[Theta]]] - Cos[py]) + 2 tdownx (Cos[k Cos[\[Theta]]] - Cos[px]) + 2 tdowny (Cos[k Sin[\[Theta]]] - Cos[py]);


		
squarePolarRadiusForRad[\[Theta]_] := Module[{\[Theta]mod},      
      \[Theta]mod=Mod[\[Theta], \[Pi]/2];
      
      If[\[Theta]mod <= \[Pi]/4,
      	\[Pi]/Cos[\[Theta]mod],
      	\[Pi]/Cos[\[Pi]/2 - \[Theta]mod]
      ]
];



HubbardTMatrixIntegral[tupx_, tdownx_, tupy_, tdowny_, px_, py_, nk_, n\[Theta]_] := Module[
	{
		\[Theta], d\[Theta], kPV, dkPV, poles, poleWeights
	},
	
	{\[Theta], d\[Theta]} = GaussianQuadratureWeights[n\[Theta], 0, 2\[Pi]]\[Transpose];
	
	(2\[Pi])^-2 Table[
	    poles = FindAllCrossings[
    		HubbardInvFreeGreenFuncPolarCoord[tupx, tdownx, tupy, tdowny, px, py, k, \[Theta]i],
    		{k, 0, Sqrt[2] \[Pi]}
    	]; (*locate all poles*)
    	
    	poles = Select[
        	poles,
        	Abs[# Cos[\[Theta]i]] < \[Pi] && Abs[# Sin[\[Theta]i]] < \[Pi] &
    	]; (*only includes poles in the integration domain -\[Pi] < kx < \[Pi], -\[Pi] < ky < \[Pi]*)
    
    	{kPV, dkPV} = setupKIntPV[nk, poles, 0, squarePolarRadiusForRad[\[Theta]i]];(*for re part*)
    
    	poleWeights=Table[
    		1 / (Abs @ ND[HubbardInvFreeGreenFuncPolarCoord[tupx, tdownx, tupy, tdowny, px, py, k, \[Theta]i], k, poles[[i]]]),
    		{i, Length @ poles}
    	]; 
	
    	dkPV . (#/(HubbardInvFreeGreenFuncPolarCoord[tupx, tdownx, tupy, tdowny, px, py, #, \[Theta]i]) &/@ kPV) - I \[Pi] Total[poles * poleWeights]

    	, {\[Theta]i, \[Theta]}
    ] . d\[Theta]

];



(* ::Section:: *)
(*4. Eigensystems*)


(* ::Text:: *)
(*	Eigensystems of the single-particle Hamiltonians*)


EigSystH0\[Sigma][v\[Sigma]_, nmax_, k_] := Module[{H0\[Sigma], esyst},
	H0\[Sigma] = DiagonalMatrix[Table[
		(k + 2. \[Pi] n)^2 / 2,
		{n, -nmax, nmax}
	] + v\[Sigma] / 2.] - (DiagonalMatrix[ConstantArray[v\[Sigma] / 4., 2 nmax], -1] + DiagonalMatrix[ConstantArray[v\[Sigma] / 4., 2 nmax], 1]);
	sortEigensystemByMagnitude @ Eigensystem @ H0\[Sigma]
];

Ek0\[Sigma][v\[Sigma]_, nmax_, k_?NumericQ] := Module[{H0\[Sigma]},
	H0\[Sigma] = DiagonalMatrix[Table[
		(k + 2. \[Pi] n)^2 / 2,
		{n, -nmax, nmax}
	] + v\[Sigma] / 2.] - (DiagonalMatrix[ConstantArray[v\[Sigma] / 4., 2 nmax], -1] + DiagonalMatrix[ConstantArray[v\[Sigma] / 4., 2 nmax], 1]);
	Eigenvalues[SparseArray[H0\[Sigma]], -1][[1]]
];

EigSysH0\[Sigma]Ground[v\[Sigma]_, nmax_, k_] := Module[{H0\[Sigma]},
	H0\[Sigma] = DiagonalMatrix[Table[
		(k + 2. \[Pi] n)^2 / 2,
		{n, -nmax, nmax}
	] + v\[Sigma] / 2.] - (DiagonalMatrix[ConstantArray[v\[Sigma] / 4., 2 nmax], -1] + DiagonalMatrix[ConstantArray[v\[Sigma] / 4., 2 nmax], 1]);
	Eigensystem[SparseArray[H0\[Sigma]], -1]
];

EigSystHc[vup_, vdown_, nmax_, k_] := Module[{Hc},
	Hc = DiagonalMatrix[Table[
		(k + 2. \[Pi] n)^2 / 4,
		{n, -nmax, nmax}
	] + (vup + vdown) / 2.] - (vup + vdown) / 4. (DiagonalMatrix[ConstantArray[1 ,2 nmax], -1] + DiagonalMatrix[ConstantArray[1, 2 nmax], 1]);
	sortEigensystemByMagnitude @ Eigensystem @ Hc
];


(* ::Section:: *)
(*5. Interacting Hamiltonian tensor*)


(* ::Text:: *)
(*	Building the interacting Hamiltonian tensor*)


buildHintTensor[vup_, vdown_, kgrid_, nCMmax_, nRelmax_] := Module[
	{
		\[Eta]eigs\[LetterSpace]\[Nu], \[Eta]ergs\[LetterSpace]\[Nu]n, \[Phi]eigs0up\[LetterSpace]\[Nu], \[Phi]ergs0up\[LetterSpace]\[Nu]n, \[Phi]eigs0up\[LetterSpace]k\[Nu],
		\[Phi]eigs0down\[LetterSpace]\[Nu], \[Phi]ergs0down\[LetterSpace]\[Nu]n, \[Nu]ergs\[LetterSpace]\[Nu]n1n2, Hintq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2,
		\[Phi]ergs0\[LetterSpace]\[Nu]upn1\[Nu]downn2, Hintq\[LetterSpace]\[Nu]\[Nu]1\[Nu]2, \[Phi]eigs0down\[LetterSpace]k\[Nu], shift\[Phi],
		\[Phi]ergsupdown\[LetterSpace]N\[Nu]1\[Nu]2, lowEnergyGridSize
	},
	
	lowEnergyGridSize = Ceiling[nCMmax/2] + nRelmax + 1;
	
	(*noninteracting dimer eigensystem*)
	{\[Eta]eigs\[LetterSpace]\[Nu], \[Eta]ergs\[LetterSpace]\[Nu]n} = Chop @ EigSystHc[vup, vdown, nCMmax, 0.]; (*factor of 2 to perform sum n_1+n_2 *)
	\[Eta]ergs\[LetterSpace]\[Nu]n = SparseArray @ \[Eta]ergs\[LetterSpace]\[Nu]n;
	
	(*while calculating matrix elements of Hintq, store the eigenvalues of the two atoms
	for later use in calculating \[CapitalPi]*)
	\[Phi]eigs0up\[LetterSpace]k\[Nu] = {}; (*actually indexed by negative k*)
	\[Phi]eigs0down\[LetterSpace]k\[Nu] = {};
	
	shift\[Phi] = (2lowEnergyGridSize + 2)/2; (*useful to perform contractions below*)


	Hintq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2 = Table[
		(*free atoms eigensystem at relative momentum k*)
		{\[Phi]eigs0up\[LetterSpace]\[Nu], \[Phi]ergs0up\[LetterSpace]\[Nu]n} = Chop @ EigSystH0\[Sigma][vup, lowEnergyGridSize, -k];
		{\[Phi]eigs0down\[LetterSpace]\[Nu], \[Phi]ergs0down\[LetterSpace]\[Nu]n} = Chop @ EigSystH0\[Sigma][vdown, lowEnergyGridSize, k];
		\[Phi]ergs0up\[LetterSpace]\[Nu]n = SparseArray @ \[Phi]ergs0up\[LetterSpace]\[Nu]n;
		\[Phi]ergs0down\[LetterSpace]\[Nu]n = SparseArray @ \[Phi]ergs0down\[LetterSpace]\[Nu]n;

		(*store eigenvalues for noninteracting atoms*)
		\[Phi]eigs0up\[LetterSpace]k\[Nu] = Append[\[Phi]eigs0up\[LetterSpace]k\[Nu], \[Phi]eigs0up\[LetterSpace]\[Nu]];
		\[Phi]eigs0down\[LetterSpace]k\[Nu] = Append[\[Phi]eigs0down\[LetterSpace]k\[Nu], \[Phi]eigs0down\[LetterSpace]\[Nu]];
		
		\[Phi]ergsupdown\[LetterSpace]N\[Nu]1\[Nu]2 = Table[
			Reverse[\[Phi]ergs0down\[LetterSpace]\[Nu]n[[;;, shift\[Phi] + Nvar/2 + If[Mod[Nvar, 2] == 1, 1/2, 0] - nRelmax ;; shift\[Phi] + Nvar/2 + If[Mod[Nvar, 2] == 1, -1/2, 0] + nRelmax]], 2] .
			\[Phi]ergs0up\[LetterSpace]\[Nu]n[[\[Nu]1, shift\[Phi] + Nvar/2 + If[Mod[Nvar, 2] == 1, +1/2, 0] - nRelmax ;; shift\[Phi] + Nvar/2 + If[Mod[Nvar, 2] == 1, -1/2, 0] + nRelmax]],
			{Nvar, -nCMmax, nCMmax}, {\[Nu]1, 1, 2*shift\[Phi] - 1}
		];

		(*Hintq\[LetterSpace]\[Nu]\[Nu]1\[Nu]2 = TensorMultiply[\[Eta]ergs\[LetterSpace]\[Nu]n, \[Phi]ergsupdown\[LetterSpace]N\[Nu]1\[Nu]2, {{2, 1}}]*)
		Hintq\[LetterSpace]\[Nu]\[Nu]1\[Nu]2 = \[Eta]ergs\[LetterSpace]\[Nu]n . \[Phi]ergsupdown\[LetterSpace]N\[Nu]1\[Nu]2,
		{k, kgrid}
	];
	
	{Hintq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2, \[Phi]eigs0up\[LetterSpace]k\[Nu], \[Phi]eigs0down\[LetterSpace]k\[Nu]}
]


(* ::Text:: *)
(*Special build of \[Nu]up = \[Nu]down = 0 Hint tensor.*)


buildHintTensor\[Nu]1\[Nu]2eq00[vup_, vdown_, k_, nCMmax_, nRelmax_, \[Eta]ergs\[LetterSpace]\[Nu]n_] := Module[
	{
		\[Phi]eigs0up\[LetterSpace]0, \[Phi]ergs0up\[LetterSpace]0n, \[Phi]eigs0down\[LetterSpace]0, \[Phi]ergs0down\[LetterSpace]0n,
		shift\[Phi], \[Phi]ergsupdown\[LetterSpace]N, lowEnergyGridSize
	},
	
	lowEnergyGridSize = Ceiling[nCMmax/2] + nRelmax + 1;
	
	{\[Phi]eigs0up\[LetterSpace]0, \[Phi]ergs0up\[LetterSpace]0n} = Chop @ EigSysH0\[Sigma]Ground[vup, lowEnergyGridSize, -k];
	{\[Phi]eigs0down\[LetterSpace]0, \[Phi]ergs0down\[LetterSpace]0n} = Chop @ EigSysH0\[Sigma]Ground[vdown, lowEnergyGridSize, k];
	shift\[Phi] = (2(Ceiling[nCMmax/2] + nRelmax + 1) + 2)/2;
	
	\[Phi]ergsupdown\[LetterSpace]N = Table[
		Reverse[\[Phi]ergs0down\[LetterSpace]0n[[1, shift\[Phi] + Nvar/2 + If[Mod[Nvar, 2] == 1, 1/2, 0] - nRelmax ;; shift\[Phi] + Nvar/2 + If[Mod[Nvar, 2] == 1, -1/2, 0] + nRelmax]]] .
		\[Phi]ergs0up\[LetterSpace]0n[[1, shift\[Phi] + Nvar/2 + If[Mod[Nvar, 2] == 1, +1/2, 0] - nRelmax ;; shift\[Phi] + Nvar/2 + If[Mod[Nvar, 2] == 1, -1/2, 0] + nRelmax]],(*check minus plus half stuff*)
		{Nvar, -nCMmax, nCMmax}
	];
	
	(*TensorMultiply[\[Eta]ergs\[LetterSpace]\[Nu]n, \[Phi]ergsupdown\[LetterSpace]N, {{2, 1}}]*)
	\[Eta]ergs\[LetterSpace]\[Nu]n . \[Phi]ergsupdown\[LetterSpace]N
	
];


(* ::Section:: *)
(*6. Pole matrix calculation*)


(* ::Text:: *)
(*	Calculating the pole component to the \[CapitalPi] matrix, arising from the \[Nu]up = \[Nu]down = 0 term.*)


(* ::Text:: *)
(*The integral is performed in a rather involved way due to difficulties of the poles in 2D. We first calculate the integral over a circle of radius \[Pi] in polar coordinates, which gives rise to \[CapitalPi]PoleContributionMatrixCircle. Then, the remainder of the integral is solved in cartesian coordinates, which gives rise to \[CapitalPi]PoleContributionMatrixSquare.*)
(**)
(*check if there is an issue with only considering poles within the circle!*)


\[CapitalPi]PoleComponent[vupx_, vdownx_, vupy_, vdowny_, pOnShellx_, pOnShelly_, nk_, n\[Theta]_, nkSquare_, nCMmax_, nRelmax_] := Module[
	{
		invFreeGreenFunc, invFreeGreenFuncPolarCoord, \[Eta]eigsx\[LetterSpace]\[Nu], \[Eta]ergsx\[LetterSpace]\[Nu]n, \[Theta], d\[Theta], \[CapitalPi]PoleContributionMatrixCircle,
		E0, poles, kPV, dkPV, Hintxq\[LetterSpace]kPV\[Nu]00, HintxqHintxq\[LetterSpace]kPV\[Nu]\[Nu]p00, Hintyq\[LetterSpace]kPV\[Nu]00, HintyqHintyq\[LetterSpace]kPV\[Nu]\[Nu]p00,
		Hintxq\[LetterSpace]kpole\[Nu]00, HintxqHintxq\[LetterSpace]kpole\[Nu]\[Nu]p00, Hintyq\[LetterSpace]kpole\[Nu]00, HintyqHintyq\[LetterSpace]kpole\[Nu]\[Nu]p00, progress,
		poleWeights, kSquare, dkSquare, kSquarePos, dkSquarePos, kySquare, dkySquare, Hintyq\[LetterSpace]k\[Nu]00,
		HintyqHintyq\[LetterSpace]k\[Nu]\[Nu]p00, HintyqHintyqdk\[LetterSpace]k\[Nu]\[Nu]p00, tempForHints, kxSquareTemp, dkxSquareTemp,
		kxSquare, dkxSquare, Hintxq\[LetterSpace]k\[Nu]00, HintxqHintxq\[LetterSpace]k\[Nu]\[Nu]p00, HintxqHintxq\[LetterSpace]kykx\[Nu]\[Nu]p00, \[Eta]eigsy\[LetterSpace]\[Nu], \[Eta]ergsy\[LetterSpace]\[Nu]n,
		integrationOverkx\[LetterSpace]ky\[Nu]\[Nu]p, \[CapitalPi]PoleContributionMatrixSquare, \[CapitalPi]PoleContributionMatrix, lowEnergyGridSize
	},
	
	lowEnergyGridSize = Ceiling[nCMmax/2] + nRelmax + 1;	
	E0 = Ek0\[Sigma][vupx, lowEnergyGridSize, -pOnShellx] + Ek0\[Sigma][vupy, lowEnergyGridSize, -pOnShelly] + Ek0\[Sigma][vdownx, lowEnergyGridSize, pOnShellx] + Ek0\[Sigma][vdowny, lowEnergyGridSize, pOnShelly]; 
	
	(*Green's functions*)
	invFreeGreenFunc[kx_, ky_] := E0 - (Ek0\[Sigma][vupx, lowEnergyGridSize, - kx] + Ek0\[Sigma][vupy, lowEnergyGridSize, -ky] + Ek0\[Sigma][vdownx, lowEnergyGridSize, kx] + Ek0\[Sigma][vdowny, lowEnergyGridSize, ky]);
	invFreeGreenFuncPolarCoord[k_, \[Theta]_] := E0 - (Ek0\[Sigma][vupx, lowEnergyGridSize, -k Cos[\[Theta]]] + Ek0\[Sigma][vupy, lowEnergyGridSize, -k Sin[\[Theta]]] + Ek0\[Sigma][vdownx, lowEnergyGridSize, k Cos[\[Theta]]] + Ek0\[Sigma][vdowny, lowEnergyGridSize, k Sin[\[Theta]]]);
	
	(*for buildHintTensor\[Nu]1\[Nu]2eq00 *)
	{\[Eta]eigsx\[LetterSpace]\[Nu], \[Eta]ergsx\[LetterSpace]\[Nu]n} = Chop @ EigSystHc[vupx, vdownx, nCMmax, 0.]; 
	\[Eta]ergsx\[LetterSpace]\[Nu]n = SparseArray @ \[Eta]ergsx\[LetterSpace]\[Nu]n;
	{\[Eta]eigsy\[LetterSpace]\[Nu], \[Eta]ergsy\[LetterSpace]\[Nu]n} = Chop @ EigSystHc[vupy, vdowny, nCMmax, 0.];
	\[Eta]ergsy\[LetterSpace]\[Nu]n = SparseArray @ \[Eta]ergsy\[LetterSpace]\[Nu]n;

	(*CIRCULAR PART*)
	
	{\[Theta], d\[Theta]} = GaussianQuadratureWeights[n\[Theta], 0, 2\[Pi]]\[Transpose];
	
	progress=0;
	
	\[CapitalPi]PoleContributionMatrixCircle = Monitor[
		ArrayFlatten[(2\[Pi])^-2 Total[Table[
			progress++;
			poles = FindAllCrossings[invFreeGreenFuncPolarCoord[k, \[Theta]i], {k, 0, \[Pi]}];
			If[Length @ poles == 0, "Error: Cannot find poles of free Green's function."]; (*either remove or make this a proper error*)
			
			(*for re part*)
			{kPV, dkPV} = setupKIntPV[nk, poles, 0, \[Pi]];
			Hintxq\[LetterSpace]kPV\[Nu]00 = buildHintTensor\[Nu]1\[Nu]2eq00[vupx, vdownx, #, nCMmax, nRelmax, \[Eta]ergsx\[LetterSpace]\[Nu]n] &/@ (kPV Cos[\[Theta]i]);
			HintxqHintxq\[LetterSpace]kPV\[Nu]\[Nu]p00 = Table[Outer[Times, Hintxq\[LetterSpace]kPV\[Nu]00[[i]], Hintxq\[LetterSpace]kPV\[Nu]00[[i]]], {i, Length @ kPV}];
			Hintyq\[LetterSpace]kPV\[Nu]00 = buildHintTensor\[Nu]1\[Nu]2eq00[vupy, vdowny, #, nCMmax, nRelmax, \[Eta]ergsy\[LetterSpace]\[Nu]n] &/@ (kPV Sin[\[Theta]i]);
			HintyqHintyq\[LetterSpace]kPV\[Nu]\[Nu]p00 = Table[Outer[Times, Hintyq\[LetterSpace]kPV\[Nu]00[[i]], Hintyq\[LetterSpace]kPV\[Nu]00[[i]]], {i, Length @ kPV}];
			
			(*for im part*)
			Hintxq\[LetterSpace]kpole\[Nu]00 = buildHintTensor\[Nu]1\[Nu]2eq00[vupx, vdownx, #, nCMmax, nRelmax, \[Eta]ergsx\[LetterSpace]\[Nu]n] &/@ (poles Cos[\[Theta]i]);
			HintxqHintxq\[LetterSpace]kpole\[Nu]\[Nu]p00 = Table[Outer[Times, Hintxq\[LetterSpace]kpole\[Nu]00[[i]], Hintxq\[LetterSpace]kpole\[Nu]00[[i]]], {i, Length @ poles}];
			Hintyq\[LetterSpace]kpole\[Nu]00 = buildHintTensor\[Nu]1\[Nu]2eq00[vupy, vdowny, #, nCMmax, nRelmax, \[Eta]ergsy\[LetterSpace]\[Nu]n] &/@ (poles Sin[\[Theta]i]);
			HintyqHintyq\[LetterSpace]kpole\[Nu]\[Nu]p00 = Table[Outer[Times, Hintyq\[LetterSpace]kpole\[Nu]00[[i]], Hintyq\[LetterSpace]kpole\[Nu]00[[i]]], {i, Length @ poles}];

			poleWeights = Table[1 / Abs @ ND[invFreeGreenFuncPolarCoord[k, \[Theta]i], k, poles[[i]]], {i, Length @ poles}];
			
			ParallelSum[TensorProduct[kPV[[i]] dkPV[[i]] HintxqHintxq\[LetterSpace]kPV\[Nu]\[Nu]p00[[i]] 1 / invFreeGreenFuncPolarCoord[kPV[[i]], \[Theta]i], HintyqHintyq\[LetterSpace]kPV\[Nu]\[Nu]p00[[i]]], {i, Length @ kPV}](*add parallel?*)
			-I \[Pi] Sum[TensorProduct[ poles[[i]] poleWeights[[i]] HintxqHintxq\[LetterSpace]kpole\[Nu]\[Nu]p00[[i]], HintyqHintyq\[LetterSpace]kpole\[Nu]\[Nu]p00[[i]]], {i, Length @ poles}],
			{\[Theta]i, \[Theta]}
		] * d\[Theta]]],
		Grid[{{Text[Style["Calculating T matrix (2/2): pole contribution ", Darker@ColorData[97, "ColorList"][[1]]]], ProgressIndicator[progress, {1, Length @ \[Theta]}]}}]
	];
	

	

	(*SQUARE PART*)

	{kSquare, dkSquare} = setupKIntReg[Round[nkSquare, 2] (*Round to ensure even integration around 0*)];
	{kSquarePos,dkSquarePos}={kSquare[[Round[nkSquare,2]/2 + 1 ;;]], dkSquare[[Round[nkSquare,2]/2 + 1 ;;]]};
	{kySquare,dkySquare}={kSquare, dkSquare};
	
	(*go through and check if kSquarePos is used*)

	Hintyq\[LetterSpace]k\[Nu]00 = buildHintTensor\[Nu]1\[Nu]2eq00[vupy, vdowny, #, nCMmax, nRelmax, \[Eta]ergsy\[LetterSpace]\[Nu]n] &/@ kySquare;
	HintyqHintyq\[LetterSpace]k\[Nu]\[Nu]p00 = Table[Outer[Times, Hintyq\[LetterSpace]k\[Nu]00[[i]], Hintyq\[LetterSpace]k\[Nu]00[[i]]], {i, Length@kySquare}];
	HintyqHintyqdk\[LetterSpace]k\[Nu]\[Nu]p00 = HintyqHintyq\[LetterSpace]k\[Nu]\[Nu]p00 * dkySquare / (2\[Pi]);

	(*only need to calculate Hints for positive ky*)
	tempForHints = ParallelTable[
		{kxSquareTemp, dkxSquareTemp} = {((\[Pi] - Sqrt[\[Pi]^2 - kyi^2]) / \[Pi] kSquarePos + Sqrt[\[Pi]^2 - kyi^2]), (\[Pi] - Sqrt[\[Pi]^2 - kyi^2]) / \[Pi] dkSquarePos};
		{kxSquare, dkxSquare} = {Join[-Reverse @ kxSquareTemp, kxSquareTemp], Join[Reverse @ dkxSquareTemp, dkxSquareTemp]};
		Hintxq\[LetterSpace]k\[Nu]00 = buildHintTensor\[Nu]1\[Nu]2eq00[vupx, vdownx, #, nCMmax, nRelmax, \[Eta]ergsx\[LetterSpace]\[Nu]n] &/@ kxSquare;
		HintxqHintxq\[LetterSpace]k\[Nu]\[Nu]p00 = Table[Outer[Times, Hintxq\[LetterSpace]k\[Nu]00[[ix]], Hintxq\[LetterSpace]k\[Nu]00[[ix]]], {ix, Length @ kxSquare}],
		{kyi, kSquarePos}];
		
	(*full Hint given by*)
	HintxqHintxq\[LetterSpace]kykx\[Nu]\[Nu]p00 = Join[Reverse @ tempForHints, tempForHints];


	integrationOverkx\[LetterSpace]ky\[Nu]\[Nu]p = Table[
		{kxSquareTemp, dkxSquareTemp} = {((\[Pi] - Sqrt[\[Pi]^2 - kyi^2]) / \[Pi] kSquarePos + Sqrt[\[Pi]^2 - kyi^2]), (\[Pi] - Sqrt[\[Pi]^2 - kyi^2]) / \[Pi] dkSquarePos};
		{kxSquare, dkxSquare} = {Join[-Reverse @ kxSquareTemp, kxSquareTemp], Join[Reverse @ dkxSquareTemp, dkxSquareTemp]};
		TensorMultiply[HintxqHintxq\[LetterSpace]kykx\[Nu]\[Nu]p00[[Position[kySquare, kyi][[1, 1]]]](1 / invFreeGreenFunc[#, kyi] &/@ kxSquare), dkxSquare / (2\[Pi]), {{1, 1}}],
		{kyi, kySquare}
	];

	\[CapitalPi]PoleContributionMatrixSquare = ArrayFlatten @ TensorMultiply[integrationOverkx\[LetterSpace]ky\[Nu]\[Nu]p, HintyqHintyqdk\[LetterSpace]k\[Nu]\[Nu]p00, {{1, 1}}];
	(*\[CapitalPi]PoleContributionMatrixSquare = ArrayFlatten @ Transpose[integrationOverkx\[LetterSpace]ky\[Nu]\[Nu]p, {3,1,2}] . HintyqHintyqdk\[LetterSpace]k\[Nu]\[Nu]p00;*)

	\[CapitalPi]PoleContributionMatrix = \[CapitalPi]PoleContributionMatrixSquare + \[CapitalPi]PoleContributionMatrixCircle
];


(* ::Section:: *)
(*7. Hubbard U calculation*)


(* ::Text:: *)
(*	Main functions of this package used to calculate the Hubbard on-site interaction energy U*)
(**)
(*Pure-2D on-shell T matrix calculation*)


(*setupTMatrixOnShellPure2D[vupx_, vdownx_, vupy_, vdowny_, nCMmax_, nRelmax_, nkReg_, nkPoleCirc_, n\[Theta]PoleCirc_, nkPoleSquare_, pOnShellx_, pOnShelly_, \[CapitalLambda]_] := Module[
	{
		E0, kReg, dkReg, Hintxq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2, Hintxdq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2, HintxqHintxdq\[LetterSpace]\[Nu]\[Nu]pk\[Nu]1\[Nu]2, \[Phi]eigs0upx\[LetterSpace]k\[Nu],
		\[Phi]eigs0downx\[LetterSpace]k\[Nu], neg\[Phi]eigs0upx\[LetterSpace]k\[Nu], neg\[Phi]eigs0downx\[LetterSpace]k\[Nu], negEigsx\[LetterSpace]k\[Nu]1\[Nu]2, negEigsPlusEnergyx\[LetterSpace]k\[Nu]1\[Nu]2,
		Hintyq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2, Hintydq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2, HintyqHintydq\[LetterSpace]\[Nu]\[Nu]pk\[Nu]1\[Nu]2, \[Phi]eigs0upy\[LetterSpace]k\[Nu], \[Phi]eigs0downy\[LetterSpace]k\[Nu], neg\[Phi]eigs0upy\[LetterSpace]k\[Nu],
		neg\[Phi]eigs0downy\[LetterSpace]k\[Nu], negEigsy\[LetterSpace]k\[Nu]1\[Nu]2, \[CapitalPi]NonPoleMatrix, kernelPoleRemoved\[LetterSpace]\[Nu]1\[Nu]2\[Nu]1\[Nu]2, \[CapitalPi]PoleMatrix,
		\[CapitalPi]Matrix, \[Eta]eigsx\[LetterSpace]\[Nu], \[Eta]ergsx\[LetterSpace]\[Nu]n, Hintxp\[LetterSpace]\[Nu]00, \[Eta]eigsy\[LetterSpace]\[Nu], \[Eta]ergsy\[LetterSpace]\[Nu]n, Hintyp\[LetterSpace]\[Nu]00, HintxpHintyp\[LetterSpace]\[Nu]\[Nu]p, 
		progress, TMatrixOnShellPure2D, lowEnergyGridSize, highEnergyIntegrals\[LetterSpace]N, \[CapitalLambda]Nx, \[CapitalLambda]Ny, highEnergyIntegrals\[LetterSpace]NxNy,
		\[Eta]ergsx\[Eta]ergsx\[LetterSpace]\[Nu]n\[LetterSpace]n\[Nu]\[Nu]p, \[Eta]ergsy\[Eta]ergsy\[LetterSpace]\[Nu]n\[LetterSpace]n\[Nu]\[Nu]p, \[CapitalPi]HighEnergyCorrection
	},
	
	
	lowEnergyGridSize = Ceiling[nCMmax/2] + nRelmax + 1;	
	E0 = Ek0\[Sigma][vupx, lowEnergyGridSize, -pOnShellx] + Ek0\[Sigma][vupy, lowEnergyGridSize, -pOnShelly] + Ek0\[Sigma][vdownx, lowEnergyGridSize, pOnShellx] + Ek0\[Sigma][vdowny, lowEnergyGridSize, pOnShelly]; 
	
	
	(*HIGH ENERGY CORRECTION*)
	{\[Eta]eigsx\[LetterSpace]\[Nu], \[Eta]ergsx\[LetterSpace]\[Nu]n} = Chop @ EigSystHc[vupx, vdownx, nCMmax, 0.]; 
	\[Eta]ergsx\[LetterSpace]\[Nu]n = SparseArray @ \[Eta]ergsx\[LetterSpace]\[Nu]n;
	\[Eta]ergsx\[Eta]ergsx\[LetterSpace]\[Nu]n\[LetterSpace]n\[Nu]\[Nu]p = Table[TensorProduct[\[Eta]ergsx\[LetterSpace]\[Nu]n[[;;, Nvar]], \[Eta]ergsx\[LetterSpace]\[Nu]n[[;;, Nvar]]], {Nvar, 1, 2nCMmax + 1}];
	{\[Eta]eigsy\[LetterSpace]\[Nu], \[Eta]ergsy\[LetterSpace]\[Nu]n} = Chop @ EigSystHc[vupy, vdowny, nCMmax, 0.]; 
	\[Eta]ergsy\[LetterSpace]\[Nu]n = SparseArray @ \[Eta]ergsy\[LetterSpace]\[Nu]n;
	\[Eta]ergsy\[Eta]ergsy\[LetterSpace]\[Nu]n\[LetterSpace]n\[Nu]\[Nu]p = Table[TensorProduct[\[Eta]ergsy\[LetterSpace]\[Nu]n[[;;, Nvar]], \[Eta]ergsy\[LetterSpace]\[Nu]n[[;;, Nvar]]], {Nvar, 1, 2nCMmax + 1}];
	
	highEnergyIntegrals\[LetterSpace]NxNy = Table[
		\[CapitalLambda]Nx = If[EvenQ @ Nvarx, \[Pi] (2nRelmax + 1), 2\[Pi] nRelmax];
		\[CapitalLambda]Ny = If[EvenQ @ Nvary, \[Pi] (2nRelmax + 1), 2\[Pi] nRelmax];
		4/(2\[Pi])^2 (NIntegrate[1/(E0 - ((2\[Pi] Nvarx )^2 + (2\[Pi] Nvary )^2)/4 - kx^2 - ky^2), {kx, \[CapitalLambda]Nx, \[CapitalLambda]}, {ky, 0, Sqrt[\[CapitalLambda]^2 - kx^2]}]
		+ NIntegrate[1/(E0 - ((2\[Pi] Nvarx )^2 + (2\[Pi] Nvary )^2)/4 - kx^2 - ky^2), {kx, 0, \[CapitalLambda]Nx}, {ky, \[CapitalLambda]Ny, Sqrt[\[CapitalLambda]^2 - kx^2]}]),
		{Nvarx, -nCMmax, nCMmax},
		{Nvary, -nCMmax, nCMmax}	(*not convinced this factor of 4 is correct for finite q*)
	];
	
	\[CapitalPi]HighEnergyCorrection = ArrayFlatten @ Activate @ TensorContract[Inactive[TensorProduct][
							\[Eta]ergsx\[Eta]ergsx\[LetterSpace]\[Nu]n\[LetterSpace]n\[Nu]\[Nu]p,
							\[Eta]ergsy\[Eta]ergsy\[LetterSpace]\[Nu]n\[LetterSpace]n\[Nu]\[Nu]p,
							highEnergyIntegrals\[LetterSpace]NxNy],
							{{1, 7}, {4, 8}}] + 1/(2\[Pi]) Log[\[CapitalLambda]/(2 \[Pi] nRelmax)]IdentityMatrix[(2nCMmax + 1)^2]; 

	
	(*NON-POLE CONTRIBUTION*)
	{kReg, dkReg} = setupKIntReg[nkReg];	
	{Hintxq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2, \[Phi]eigs0upx\[LetterSpace]k\[Nu], \[Phi]eigs0downx\[LetterSpace]k\[Nu]} = buildHintTensor[vupx, vdownx, kReg, nCMmax, nRelmax];
	Hintxdq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2 = dkReg / (2\[Pi]) Hintxq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2;
	HintxqHintxdq\[LetterSpace]\[Nu]\[Nu]pk\[Nu]1\[Nu]2 = Table[Hintxq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2[[;;, \[Nu], ;;]] Hintxdq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2[[;;, \[Nu]p, ;;]], {\[Nu], 2nCMmax + 1}, {\[Nu]p, 2nCMmax + 1}];
	neg\[Phi]eigs0upx\[LetterSpace]k\[Nu] = -\[Phi]eigs0upx\[LetterSpace]k\[Nu];
	neg\[Phi]eigs0downx\[LetterSpace]k\[Nu] = -\[Phi]eigs0downx\[LetterSpace]k\[Nu];
	negEigsx\[LetterSpace]k\[Nu]1\[Nu]2 = Table[Outer[(#1 + #2) &, neg\[Phi]eigs0upx\[LetterSpace]k\[Nu][[kxi]], neg\[Phi]eigs0downx\[LetterSpace]k\[Nu][[kxi]]], {kxi, Length @ kReg}];
	negEigsPlusEnergyx\[LetterSpace]k\[Nu]1\[Nu]2 = negEigsx\[LetterSpace]k\[Nu]1\[Nu]2 + E0;
	
	{Hintyq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2, \[Phi]eigs0upy\[LetterSpace]k\[Nu], \[Phi]eigs0downy\[LetterSpace]k\[Nu]} = buildHintTensor[vupy, vdowny, kReg, nCMmax, nRelmax];
	Hintydq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2 = dkReg / (2\[Pi]) Hintyq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2;
	HintyqHintydq\[LetterSpace]\[Nu]\[Nu]pk\[Nu]1\[Nu]2 = Table[Hintyq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2[[;;, \[Nu], ;;]] Hintydq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2[[;;, \[Nu]p, ;;]], {\[Nu], 2nCMmax + 1}, {\[Nu]p, 2nCMmax + 1}];
	neg\[Phi]eigs0upy\[LetterSpace]k\[Nu] = -\[Phi]eigs0upy\[LetterSpace]k\[Nu];
	neg\[Phi]eigs0downy\[LetterSpace]k\[Nu] = -\[Phi]eigs0downy\[LetterSpace]k\[Nu];
	negEigsy\[LetterSpace]k\[Nu]1\[Nu]2 = Table[Outer[(#1 + #2) &, neg\[Phi]eigs0upy\[LetterSpace]k\[Nu][[kyi]], neg\[Phi]eigs0downy\[LetterSpace]k\[Nu][[kyi]]], {kyi, Length @ kReg}];
	
	SetSharedVariable[progress];
	progress=0;
	\[CapitalPi]NonPoleMatrix = Monitor[
		ArrayFlatten @ ParallelSum[progress++;
			kernelPoleRemoved\[LetterSpace]\[Nu]1\[Nu]2\[Nu]1\[Nu]2 = Outer[1 / (#1 + #2) &, negEigsPlusEnergyx\[LetterSpace]k\[Nu]1\[Nu]2[[kxi]], negEigsy\[LetterSpace]k\[Nu]1\[Nu]2[[kyi]]];
			kernelPoleRemoved\[LetterSpace]\[Nu]1\[Nu]2\[Nu]1\[Nu]2[[1,1,1,1]] = 0; (*remove divergence*)
			Activate @ TensorContract[Inactive[TensorProduct][
				HintxqHintxdq\[LetterSpace]\[Nu]\[Nu]pk\[Nu]1\[Nu]2[[;;, ;;, kxi, ;;, ;;]],
				HintyqHintydq\[LetterSpace]\[Nu]\[Nu]pk\[Nu]1\[Nu]2[[;;, ;;, kyi, ;;, ;;]],
				kernelPoleRemoved\[LetterSpace]\[Nu]1\[Nu]2\[Nu]1\[Nu]2],
				{{3, 9}, {4, 10}, {7, 11}, {8, 12}}],
			{kxi, Length @ kReg},
			{kyi, Length @ kReg}
		] + \[CapitalPi]HighEnergyCorrection (*note this addition!!*),
		Grid[{{Text[Style["Calculating T matrix (1/2): regular contribution ", Darker@ColorData[97, "ColorList"][[1]]]], ProgressIndicator[progress, {1, (Length @ kReg)^2}]}}]
	];
	
	
	(*POLE CONTRIBUTION*)
	\[CapitalPi]PoleMatrix = \[CapitalPi]PoleComponent[vupx, vdownx, vupy, vdowny, pOnShellx, pOnShelly, nkPoleCirc, n\[Theta]PoleCirc, nkPoleSquare, nCMmax, nRelmax];
	
	(*T Matrix*)
	\[CapitalPi]Matrix = \[CapitalPi]NonPoleMatrix + \[CapitalPi]PoleMatrix;
		

	Hintxp\[LetterSpace]\[Nu]00 = buildHintTensor\[Nu]1\[Nu]2eq00[vupx, vdownx, pOnShellx, nCMmax, nRelmax, \[Eta]ergsx\[LetterSpace]\[Nu]n];
	Hintyp\[LetterSpace]\[Nu]00 = buildHintTensor\[Nu]1\[Nu]2eq00[vupy, vdowny, pOnShelly, nCMmax, nRelmax, \[Eta]ergsy\[LetterSpace]\[Nu]n];	
	
	HintxpHintyp\[LetterSpace]\[Nu]\[Nu]p = Flatten @ Outer[Times, Hintxp\[LetterSpace]\[Nu]00, Hintyp\[LetterSpace]\[Nu]00];

	TMatrixOnShellPure2D[logadinv_] := HintxpHintyp\[LetterSpace]\[Nu]\[Nu]p . Inverse[-1 / (2\[Pi]) (-logadinv + Log[2\[Pi] nRelmax]) IdentityMatrix[Length @ \[CapitalPi]Matrix] - \[CapitalPi]Matrix] . HintxpHintyp\[LetterSpace]\[Nu]\[Nu]p;
	TMatrixOnShellPure2D
]*)


(* ::Text:: *)
(*Quasi-2D on-shell T matrix calculation*)


setupTMatrixOnShellQuasi2D[vupx_, vdownx_, vupy_, vdowny_, nCMmax_, nRelmax_, nkReg_, nkPoleCirc_, n\[Theta]PoleCirc_, nkPoleSquare_, pOnShellx_, pOnShelly_, \[CapitalLambda]_, \[Omega]z_] := Module[
	{
		E0, kReg, dkReg, Hintxq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2, Hintxdq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2, HintxqHintxdq\[LetterSpace]\[Nu]\[Nu]pk\[Nu]1\[Nu]2, \[Phi]eigs0upx\[LetterSpace]k\[Nu],
		\[Phi]eigs0downx\[LetterSpace]k\[Nu], neg\[Phi]eigs0upx\[LetterSpace]k\[Nu], neg\[Phi]eigs0downx\[LetterSpace]k\[Nu], negEigsx\[LetterSpace]k\[Nu]1\[Nu]2, negEigsPlusEnergyx\[LetterSpace]k\[Nu]1\[Nu]2,
		Hintyq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2, Hintydq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2, HintyqHintydq\[LetterSpace]\[Nu]\[Nu]pk\[Nu]1\[Nu]2, \[Phi]eigs0upy\[LetterSpace]k\[Nu], \[Phi]eigs0downy\[LetterSpace]k\[Nu], neg\[Phi]eigs0upy\[LetterSpace]k\[Nu],
		neg\[Phi]eigs0downy\[LetterSpace]k\[Nu], negEigsy\[LetterSpace]k\[Nu]1\[Nu]2, \[CapitalPi]NonPoleMatrix, kernelPoleRemoved\[LetterSpace]\[Nu]1\[Nu]2\[Nu]1\[Nu]2, \[CapitalPi]PoleMatrix,
		\[CapitalPi]Matrix, \[Eta]eigsx\[LetterSpace]\[Nu], \[Eta]ergsx\[LetterSpace]\[Nu]n, Hintxp\[LetterSpace]\[Nu]00, \[Eta]eigsy\[LetterSpace]\[Nu], \[Eta]ergsy\[LetterSpace]\[Nu]n, Hintyp\[LetterSpace]\[Nu]00, HintxpHintyp\[LetterSpace]\[Nu]\[Nu]p, 
		progress, TMatrixOnShellQuasi2D, counterTerm, counterTermMatrix, lowEnergyGridSize,
		\[CapitalPi]HighEnergyCorrection, \[Eta]ergsx\[Eta]ergsx\[LetterSpace]\[Nu]n\[LetterSpace]n\[Nu]\[Nu]p, \[Eta]ergsy\[Eta]ergsy\[LetterSpace]\[Nu]n\[LetterSpace]n\[Nu]\[Nu]p, highEnergyIntegrals\[LetterSpace]NxNy,
		\[CapitalLambda]Nx, \[CapitalLambda]Ny, smallestPotPos, lowestMomentumPos, combinedEigs, negEigsPlusEnergy\[LetterSpace]\[Nu]1\[Nu]2\[Nu]1\[Nu]2, 
		quasi2DSumApproxPos, \[CapitalLambda]1
	},
	
	
	lowEnergyGridSize = Ceiling[nCMmax/2] + nRelmax + 1;	
	E0 = Ek0\[Sigma][vupx, lowEnergyGridSize, -pOnShellx] + Ek0\[Sigma][vupy, lowEnergyGridSize, -pOnShelly] + Ek0\[Sigma][vdownx, lowEnergyGridSize, pOnShellx] + Ek0\[Sigma][vdowny, lowEnergyGridSize, pOnShelly]; 
	
	{\[Eta]eigsx\[LetterSpace]\[Nu], \[Eta]ergsx\[LetterSpace]\[Nu]n} = Chop @ EigSystHc[vupx, vdownx, nCMmax, 0.]; 
	\[Eta]ergsx\[LetterSpace]\[Nu]n = SparseArray @ \[Eta]ergsx\[LetterSpace]\[Nu]n;
	\[Eta]ergsx\[Eta]ergsx\[LetterSpace]\[Nu]n\[LetterSpace]n\[Nu]\[Nu]p = Table[TensorProduct[\[Eta]ergsx\[LetterSpace]\[Nu]n[[;;, Nvar]], \[Eta]ergsx\[LetterSpace]\[Nu]n[[;;, Nvar]]], {Nvar, 1, 2nCMmax + 1}];
	{\[Eta]eigsy\[LetterSpace]\[Nu], \[Eta]ergsy\[LetterSpace]\[Nu]n} = Chop @ EigSystHc[vupy, vdowny, nCMmax, 0.]; 
	\[Eta]ergsy\[LetterSpace]\[Nu]n = SparseArray @ \[Eta]ergsy\[LetterSpace]\[Nu]n;
	\[Eta]ergsy\[Eta]ergsy\[LetterSpace]\[Nu]n\[LetterSpace]n\[Nu]\[Nu]p = Table[TensorProduct[\[Eta]ergsy\[LetterSpace]\[Nu]n[[;;, Nvar]], \[Eta]ergsy\[LetterSpace]\[Nu]n[[;;, Nvar]]], {Nvar, 1, 2nCMmax + 1}];
	
	\[CapitalLambda]1=Sqrt[2]\[Pi] (2nRelmax + 1);
	highEnergyIntegrals\[LetterSpace]NxNy = Table[
		\[CapitalLambda]Nx = If[EvenQ @ Nvarx, \[Pi] (2nRelmax + 1), 2\[Pi] nRelmax];
		\[CapitalLambda]Ny = If[EvenQ @ Nvary, \[Pi] (2nRelmax + 1), 2\[Pi] nRelmax];
		4/(2\[Pi])^2 (NIntegrate[quasi2DSum0toInf[E0 - ((2\[Pi] Nvarx )^2 + (2\[Pi] Nvary )^2)/4 - kx^2 - ky^2, -\[Omega]z], {kx, \[CapitalLambda]Nx, \[CapitalLambda]1}, {ky, 0, Sqrt[\[CapitalLambda]1^2 - kx^2]}]
		+ NIntegrate[quasi2DSum0toInf[E0 - ((2\[Pi] Nvarx )^2 + (2\[Pi] Nvary )^2)/4 - kx^2 - ky^2, -\[Omega]z], {kx, 0, \[CapitalLambda]Nx}, {ky, \[CapitalLambda]Ny, Sqrt[\[CapitalLambda]1^2 - kx^2]}])
		+1/(2\[Pi]) NIntegrate[k quasi2DSum0toInf[E0 - ((2\[Pi] Nvarx )^2 + (2\[Pi] Nvary )^2)/4 - k^2, -\[Omega]z],{k,\[CapitalLambda]1,\[CapitalLambda]}],
		{Nvarx, -nCMmax, nCMmax},
		{Nvary, -nCMmax, nCMmax}	
	];  (*check this high energy correction at some point for any mistakes*)
	

	
	\[CapitalPi]HighEnergyCorrection = ArrayFlatten @ Activate @ TensorContract[Inactive[TensorProduct][
							\[Eta]ergsx\[Eta]ergsx\[LetterSpace]\[Nu]n\[LetterSpace]n\[Nu]\[Nu]p,
							\[Eta]ergsy\[Eta]ergsy\[LetterSpace]\[Nu]n\[LetterSpace]n\[Nu]\[Nu]p,
							highEnergyIntegrals\[LetterSpace]NxNy],
							{{1, 7}, {4, 8}}] + 1/(2\[Pi]) NIntegrate[k quasi2DSum0toInf[k^2, \[Omega]z], {k, 2 \[Pi] nRelmax, \[CapitalLambda]}] IdentityMatrix[(2nCMmax + 1)^2]; 
	
	counterTermMatrix = 1/(2\[Pi]) NIntegrate[k quasi2DSum1toInf[k^2, \[Omega]z], {k, 0, 2\[Pi] nRelmax}] IdentityMatrix[(2nCMmax + 1)^2];	
	
	(*NON-POLE CONTRIBUTION*)
	{kReg, dkReg} = setupKIntReg[nkReg];	
	{Hintxq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2, \[Phi]eigs0upx\[LetterSpace]k\[Nu], \[Phi]eigs0downx\[LetterSpace]k\[Nu]} = buildHintTensor[vupx, vdownx, kReg, nCMmax, nRelmax];
	Hintxdq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2 = dkReg / (2\[Pi]) Hintxq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2;
	HintxqHintxdq\[LetterSpace]\[Nu]\[Nu]pk\[Nu]1\[Nu]2 = Table[Hintxq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2[[;;, \[Nu], ;;]] Hintxdq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2[[;;, \[Nu]p, ;;]], {\[Nu], 2nCMmax + 1}, {\[Nu]p, 2nCMmax + 1}];
	neg\[Phi]eigs0upx\[LetterSpace]k\[Nu] = -\[Phi]eigs0upx\[LetterSpace]k\[Nu];
	neg\[Phi]eigs0downx\[LetterSpace]k\[Nu] = -\[Phi]eigs0downx\[LetterSpace]k\[Nu];
	negEigsx\[LetterSpace]k\[Nu]1\[Nu]2 = Table[Outer[(#1 + #2) &, neg\[Phi]eigs0upx\[LetterSpace]k\[Nu][[kxi]], neg\[Phi]eigs0downx\[LetterSpace]k\[Nu][[kxi]]], {kxi, Length @ kReg}];
	negEigsPlusEnergyx\[LetterSpace]k\[Nu]1\[Nu]2 = negEigsx\[LetterSpace]k\[Nu]1\[Nu]2 + E0;
	
	{Hintyq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2, \[Phi]eigs0upy\[LetterSpace]k\[Nu], \[Phi]eigs0downy\[LetterSpace]k\[Nu]} = buildHintTensor[vupy, vdowny, kReg, nCMmax, nRelmax];
	Hintydq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2 = dkReg / (2\[Pi]) Hintyq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2;
	HintyqHintydq\[LetterSpace]\[Nu]\[Nu]pk\[Nu]1\[Nu]2 = Table[Hintyq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2[[;;, \[Nu], ;;]] Hintydq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2[[;;, \[Nu]p, ;;]], {\[Nu], 2nCMmax + 1}, {\[Nu]p, 2nCMmax + 1}];
	neg\[Phi]eigs0upy\[LetterSpace]k\[Nu] = -\[Phi]eigs0upy\[LetterSpace]k\[Nu];
	neg\[Phi]eigs0downy\[LetterSpace]k\[Nu] = -\[Phi]eigs0downy\[LetterSpace]k\[Nu];
	negEigsy\[LetterSpace]k\[Nu]1\[Nu]2 = Table[Outer[(#1 + #2) &, neg\[Phi]eigs0upy\[LetterSpace]k\[Nu][[kyi]], neg\[Phi]eigs0downy\[LetterSpace]k\[Nu][[kyi]]], {kyi, Length @ kReg}];
	
	
	(*find position to begin applying approximation to quasi-2D sum*)
	smallestPotPos = FirstPosition[{vupx, vdownx, vupy, vdowny}, Min @ {vupx, vdownx, vupy, vdowny}][[1]];
	lowestMomentumPos = Ceiling @ ((nkReg+1)/2);
	combinedEigs = {\[Phi]eigs0upx\[LetterSpace]k\[Nu], \[Phi]eigs0downx\[LetterSpace]k\[Nu], \[Phi]eigs0upy\[LetterSpace]k\[Nu], \[Phi]eigs0downy\[LetterSpace]k\[Nu]};
	quasi2DSumApproxPos = FirstPosition[1/(2 \[Omega]z) (combinedEigs[[smallestPotPos, lowestMomentumPos, ;;]] + Total@combinedEigs[[DeleteCases[Range[1, 4], x_ /; x == smallestPotPos], lowestMomentumPos, 1]]), x_ /; x > 4, {-1}][[1]]; (*ignoring energy AND CHANGED FROM 4*)
	If[EvenQ[quasi2DSumApproxPos], quasi2DSumApproxPos++];
	
	
	SetSharedVariable[progress];
	progress=0;
	\[CapitalPi]NonPoleMatrix = Monitor[
		ArrayFlatten @ ParallelSum[
			progress++;
			negEigsPlusEnergy\[LetterSpace]\[Nu]1\[Nu]2\[Nu]1\[Nu]2 = Outer[Plus, negEigsPlusEnergyx\[LetterSpace]k\[Nu]1\[Nu]2[[kxi]], negEigsy\[LetterSpace]k\[Nu]1\[Nu]2[[kyi]]];
			
			If[quasi2DSumApproxPos == -1, 
				kernelPoleRemoved\[LetterSpace]\[Nu]1\[Nu]2\[Nu]1\[Nu]2 = quasi2DSum0toInf[negEigsPlusEnergy\[LetterSpace]\[Nu]1\[Nu]2\[Nu]1\[Nu]2, -\[Omega]z];,
				kernelPoleRemoved\[LetterSpace]\[Nu]1\[Nu]2\[Nu]1\[Nu]2 = quasi2DSum0toInfApprox[negEigsPlusEnergy\[LetterSpace]\[Nu]1\[Nu]2\[Nu]1\[Nu]2, -\[Omega]z];	
				kernelPoleRemoved\[LetterSpace]\[Nu]1\[Nu]2\[Nu]1\[Nu]2[[;; quasi2DSumApproxPos, ;; quasi2DSumApproxPos, ;; quasi2DSumApproxPos, ;; quasi2DSumApproxPos]] = 
				quasi2DSum0toInf[negEigsPlusEnergy\[LetterSpace]\[Nu]1\[Nu]2\[Nu]1\[Nu]2[[;; quasi2DSumApproxPos, ;; quasi2DSumApproxPos, ;; quasi2DSumApproxPos, ;; quasi2DSumApproxPos]], -\[Omega]z];
			];
			
			kernelPoleRemoved\[LetterSpace]\[Nu]1\[Nu]2\[Nu]1\[Nu]2[[1, 1, 1, 1]] = quasi2DSum1toInf[negEigsPlusEnergyx\[LetterSpace]k\[Nu]1\[Nu]2[[kxi, 1, 1]] + negEigsy\[LetterSpace]k\[Nu]1\[Nu]2[[kyi, 1, 1]], -\[Omega]z];
			Activate @ TensorContract[Inactive[TensorProduct][
				HintxqHintxdq\[LetterSpace]\[Nu]\[Nu]pk\[Nu]1\[Nu]2[[;;, ;;, kxi, ;;, ;;]],
				HintyqHintydq\[LetterSpace]\[Nu]\[Nu]pk\[Nu]1\[Nu]2[[;;, ;;, kyi, ;;, ;;]],
				kernelPoleRemoved\[LetterSpace]\[Nu]1\[Nu]2\[Nu]1\[Nu]2],
				{{3, 9}, {4, 10}, {7, 11}, {8, 12}}],
			{kxi, Length @ kReg},
			{kyi, Length @ kReg}
		],
		Grid[{{Text[Style["Calculating T matrix (1/2): regular contribution ", Darker@ColorData[97, "ColorList"][[1]]]], ProgressIndicator[progress, {1, (Length @ kReg)^2}]}}]
	];
	

	(*POLE CONTRIBUTION*)
	\[CapitalPi]PoleMatrix = \[CapitalPi]PoleComponent[vupx, vdownx, vupy, vdowny, pOnShellx, pOnShelly, nkPoleCirc, n\[Theta]PoleCirc, nkPoleSquare, nCMmax, nRelmax];

	(*T Matrix*)
	\[CapitalPi]Matrix = \[CapitalPi]NonPoleMatrix + \[CapitalPi]PoleMatrix + \[CapitalPi]HighEnergyCorrection + counterTermMatrix;

	Hintxp\[LetterSpace]\[Nu]00 = buildHintTensor\[Nu]1\[Nu]2eq00[vupx, vdownx, pOnShellx, nCMmax, nRelmax, \[Eta]ergsx\[LetterSpace]\[Nu]n];
	Hintyp\[LetterSpace]\[Nu]00 = buildHintTensor\[Nu]1\[Nu]2eq00[vupy, vdowny, pOnShelly, nCMmax, nRelmax, \[Eta]ergsy\[LetterSpace]\[Nu]n];	
	HintxpHintyp\[LetterSpace]\[Nu]\[Nu]p = Flatten @ Outer[Times, Hintxp\[LetterSpace]\[Nu]00, Hintyp\[LetterSpace]\[Nu]00];

	TMatrixOnShellQuasi2D[logadinv_] := HintxpHintyp\[LetterSpace]\[Nu]\[Nu]p . Inverse[-1 / (2\[Pi]) (-logadinv + Log[2 \[Pi] nRelmax]) IdentityMatrix[Length @ \[CapitalPi]Matrix] - \[CapitalPi]Matrix] . HintxpHintyp\[LetterSpace]\[Nu]\[Nu]p; 
	TMatrixOnShellQuasi2D
]


(* ::Text:: *)
(*Pure 2D Hubbard U calculation*)


(*setupHubbardUPure2D[vupx_, vdownx_, vupy_, vdowny_, nCMmax_, nRelmax_, nkReg_, nkPoleCirc_, n\[Theta]PoleCirc_, nkPoleSquare_, pOnShellx_, pOnShelly_, nkHubbard_, n\[Theta]Hubbard_, \[CapitalLambda]_] := Module[
	{
		tupx, tdownx, tupy, tdowny, hubbardTMatrixInt, effMassHubbard, TMatrixOnShellPure2D, hubbardUPure2D,
		lowEnergyGridSize
	},
	
	lowEnergyGridSize = Ceiling[nCMmax/2] + nRelmax + 1;	
	
	TMatrixOnShellPure2D = setupTMatrixOnShellPure2D[vupx, vdownx, vupy, vdowny, nCMmax, nRelmax, nkReg, nkPoleCirc, n\[Theta]PoleCirc, nkPoleSquare, pOnShellx, pOnShelly, \[CapitalLambda]];
	
	tupx = -1 / (2\[Pi]) NIntegrate[Ek0\[Sigma][vupx, lowEnergyGridSize, k] Cos[k], {k, -\[Pi], \[Pi]}];
	tdownx = -1 / (2\[Pi]) NIntegrate[Ek0\[Sigma][vdownx, lowEnergyGridSize, k] Cos[k], {k, -\[Pi], \[Pi]}];
	tupy = -1 / (2\[Pi]) NIntegrate[Ek0\[Sigma][vupy, lowEnergyGridSize, k] Cos[k], {k, -\[Pi], \[Pi]}];
	tdowny = -1 / (2\[Pi]) NIntegrate[Ek0\[Sigma][vdowny, lowEnergyGridSize, k] Cos[k], {k, -\[Pi], \[Pi]}];
	
	hubbardTMatrixInt = HubbardTMatrixIntegral[tupx, tdownx, tupy, tdowny, pOnShellx, pOnShelly, nkHubbard, n\[Theta]Hubbard];
	
	effMassHubbard = Im[-hubbardTMatrixInt] / Im[TMatrixOnShellPure2D[0]^-1];
	
	hubbardUPure2D[logadinv_] := (effMassHubbard Re[TMatrixOnShellPure2D[logadinv]^-1] + Re[hubbardTMatrixInt])^-1;(*definitely plus here for Re Hubbard? I think it is correct...*)
	
	hubbardUPure2D
	
]*)


(* ::Text:: *)
(*Quasi-2D Hubbard U calculation*)


setupHubbardUQuasi2D[vupx_, vdownx_, vupy_, vdowny_, nCMmax_, nRelmax_, nkReg_, nkPoleCirc_, n\[Theta]PoleCirc_, nkPoleSquare_, pOnShellx_, pOnShelly_, \[Omega]z_, nkHubbard_, n\[Theta]Hubbard_, \[CapitalLambda]_] := Module[
	{
		tupx, tdownx, tupy, tdowny, hubbardTMatrixInt, effMassHubbard, TMatrixOnShellQuasi2D, hubbardUQuasi2D,
		lowEnergyGridSize
	},
	
	lowEnergyGridSize = Ceiling[nCMmax/2] + nRelmax + 1;	
	
	TMatrixOnShellQuasi2D = setupTMatrixOnShellQuasi2D[vupx, vdownx, vupy, vdowny, nCMmax, nRelmax, nkReg, nkPoleCirc, n\[Theta]PoleCirc, nkPoleSquare, pOnShellx, pOnShelly, \[CapitalLambda], \[Omega]z];
	
	tupx = -1 / (2\[Pi]) NIntegrate[Ek0\[Sigma][vupx, lowEnergyGridSize, k] Cos[k], {k, -\[Pi], \[Pi]}];
	tdownx = -1 / (2\[Pi]) NIntegrate[Ek0\[Sigma][vdownx, lowEnergyGridSize, k] Cos[k], {k, -\[Pi], \[Pi]}];
	tupy = -1 / (2\[Pi]) NIntegrate[Ek0\[Sigma][vupy, lowEnergyGridSize, k] Cos[k], {k, -\[Pi], \[Pi]}];
	tdowny = -1 / (2\[Pi]) NIntegrate[Ek0\[Sigma][vdowny, lowEnergyGridSize, k] Cos[k], {k, -\[Pi], \[Pi]}];
	
	hubbardTMatrixInt = HubbardTMatrixIntegral[tupx, tdownx, tupy, tdowny, pOnShellx, pOnShelly, nkHubbard, n\[Theta]Hubbard];
	
	effMassHubbard = Im[-hubbardTMatrixInt] / Im[TMatrixOnShellQuasi2D[0]^-1];
	
	hubbardUQuasi2D[logadinv_] := (effMassHubbard Re[TMatrixOnShellQuasi2D[logadinv]^-1] + Re[hubbardTMatrixInt])^-1;
	
	hubbardUQuasi2D
]


(* ::Text:: *)
(*Quasi-2D Hubbard U calculation with automated convergence*)


setupHubbardUQuasi2D[vupx_, vdownx_, vupy_, vdowny_, pOnShellx_, pOnShelly_, \[Omega]z_, convParam_] := setupHubbardUQuasi2D[vupx, vdownx, vupy, vdowny, convParam * 2 + 3, convParam * 3 + 5, convParam * 3 + 6, 8 + convParam * 5, 8 + convParam * 5, 8 + convParam * 5, pOnShellx, pOnShelly, \[Omega]z,  50 + convParam * 10, 50 + convParam * 10, 1000 + convParam * 1000]
