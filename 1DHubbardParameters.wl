(* ::Package:: *)

(* ::Title:: *)
(*1DHubbardParameters Package*)


(* ::Subtitle:: *)
(*A package for determining the Hubbard onsite interaction energy U for a quasi-1D optical lattice*)


(* ::Author:: *)
(*Haydn S. Adlong*)


(* ::Section:: *)
(*Introduction*)


(*temp tempte mp *)
(*BeginPackage[ "1DHubbardParameters`" ]*)

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

meshgrid[x_List, y_List] := {ConstantArray[x, Length[y]], Transpose @ ConstantArray[y, Length[x]]};

TensorMultiply[A_, B_, pairs_] := Activate @ TensorContract[Inactive[TensorProduct][A, B], (# + {0, TensorRank[A]}) &/@ pairs];
(*See Brent answer: https://mathematica.stackexchange.com/questions/17758/ways-to-compute-inner-products-of-tensors*)

intervalInverse[Interval[int___]] := Interval @@ Partition[Flatten @ {int} /. {{-\[Infinity], mid___, \[Infinity]} :> {mid}, {-\[Infinity], mid__} :> {mid, \[Infinity]}, {mid__, \[Infinity]} :> {-\[Infinity], mid}, {mid___} :> {-\[Infinity], mid, \[Infinity]}}, 2]
intervalComplement[a_Interval, b__Interval] := IntervalIntersection[a, intervalInverse @ IntervalUnion[b]];
(*See Szabolcs answer: https://mathematica.stackexchange.com/questions/11345/can-mathematica-handle-open-intervals-interval-complements*)

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
(*See Stan Wagon's "Mathematica in Action" and J. M.'s persistent exhaustion answer: https://mathematica.stackexchange.com/questions/5663/about-multi-root-search-in-mathematica-for-transcendental-equations*)

quasi1DSum[a_, b_, smin_, smaxI_] := (PolyGamma[0, 1 + a/(2 b) + smaxI] - PolyGamma[0, a/(2 b) + smin])/(2 b); (*Sum[1/(a+b 2 s),{s,smin,smax}] = (PolyGamma[0,1+a/(2 b)+smax]-PolyGamma[0,a/(2 b)+smin])/(2 b)*)



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
(*	Here we perform the integral required to calculated the Hubbard model T matrix. *)


HubbardInvFreeGreenFunc[tup_, tdown_, p_, k_] :=  2 tup (Cos[-k] - Cos[-p]) + 2 tdown(Cos[k] - Cos[p])

HubbardTMatrixIntegral[tup_, tdown_, p_] := Module[ 
	{
		poles, poleWeights
	},
		
   poles = FindAllCrossings[
    		HubbardInvFreeGreenFunc[tup, tdown, p, k],
    		{k, -\[Pi],  \[Pi]}
	];

	I/2 Total[1/Abs[ND[HubbardInvFreeGreenFunc[tup, tdown, p, k], k, #]] &/@ poles]


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

		Hintq\[LetterSpace]\[Nu]\[Nu]1\[Nu]2 = TensorMultiply[\[Eta]ergs\[LetterSpace]\[Nu]n, \[Phi]ergsupdown\[LetterSpace]N\[Nu]1\[Nu]2, {{2, 1}}],
		{k, kgrid}
	];
	
	{Hintq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2, \[Phi]eigs0up\[LetterSpace]k\[Nu], \[Phi]eigs0down\[LetterSpace]k\[Nu]}
]


(* ::Text:: *)
(*Special build of \[Nu]up = \[Nu]down = 0 Hint tensor.*)


buildHintTensor\[Nu]1\[Nu]2eq00[vup_, vdown_, k_?NumericQ, nCMmax_, nRelmax_, \[Eta]ergs\[LetterSpace]\[Nu]n_] := Module[
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
	
	TensorMultiply[\[Eta]ergs\[LetterSpace]\[Nu]n, \[Phi]ergsupdown\[LetterSpace]N, {{2, 1}}]
	
];


(* ::Section:: *)
(*6. Pole matrix calculation*)


(* ::Text:: *)
(*	Calculating the pole component to the \[CapitalPi] matrix, arising from the \[Nu]up = \[Nu]down = 0 term.*)


\[CapitalPi]PoleComponent[vup_, vdown_, pOnShell_, nk_, nCMmax_, nRelmax_] := Module[
	{
		E0, invFreeGreenFunc, \[Eta]eigs\[LetterSpace]\[Nu], \[Eta]ergs\[LetterSpace]\[Nu]n, Hintq\[LetterSpace]kPV\[Nu]00, HintqHintq\[LetterSpace]kPV\[Nu]\[Nu]p00,
		Hintq\[LetterSpace]kpole\[Nu]00, HintqHintq\[LetterSpace]kpole\[Nu]\[Nu]p00, poleWeights, lowEnergyGridSize, poles,
		kPV, dkPV
	},
	
	lowEnergyGridSize = Ceiling[nCMmax/2] + nRelmax + 1;
	
	E0 = Ek0\[Sigma][vup, lowEnergyGridSize, - pOnShell] + Ek0\[Sigma][vdown, lowEnergyGridSize, pOnShell]; 
	
	(*Green's functions*)
	invFreeGreenFunc[k_] := E0 - (Ek0\[Sigma][vup, lowEnergyGridSize, - k] + Ek0\[Sigma][vdown, lowEnergyGridSize, k]);
	
	(*for buildHintTensor\[Nu]1\[Nu]2eq00 *)
	{\[Eta]eigs\[LetterSpace]\[Nu], \[Eta]ergs\[LetterSpace]\[Nu]n} = Chop @ EigSystHc[vup, vdown, nCMmax, 0];
	\[Eta]ergs\[LetterSpace]\[Nu]n = SparseArray @ \[Eta]ergs\[LetterSpace]\[Nu]n;
	
	poles = FindAllCrossings[invFreeGreenFunc[k], {k, -\[Pi], \[Pi]}];
	
	Total[{
		(*re part*)
		{kPV, dkPV} = setupKIntPV[nk, poles, -\[Pi], \[Pi]];
		Hintq\[LetterSpace]kPV\[Nu]00 = buildHintTensor\[Nu]1\[Nu]2eq00[vup, vdown, #, nCMmax, nRelmax, \[Eta]ergs\[LetterSpace]\[Nu]n] &/@ (kPV );
		HintqHintq\[LetterSpace]kPV\[Nu]\[Nu]p00 = Table[Outer[Times, Hintq\[LetterSpace]kPV\[Nu]00[[i]], Hintq\[LetterSpace]kPV\[Nu]00[[i]]], {i, Length @ kPV}];
		TensorMultiply[HintqHintq\[LetterSpace]kPV\[Nu]\[Nu]p00 (1 / invFreeGreenFunc[#] &/@ kPV), dkPV/(2\[Pi]), {{1, 1}}],
		
		(*im part*)
		Hintq\[LetterSpace]kpole\[Nu]00 = buildHintTensor\[Nu]1\[Nu]2eq00[vup, vdown, #, nCMmax, nRelmax, \[Eta]ergs\[LetterSpace]\[Nu]n] &/@ poles; (*removed \[Nu]p index*)
		HintqHintq\[LetterSpace]kpole\[Nu]\[Nu]p00 = Table[Outer[Times, Hintq\[LetterSpace]kpole\[Nu]00[[i]], Hintq\[LetterSpace]kpole\[Nu]00[[i]]], {i, Length @ poles}];
		poleWeights = (1 / Abs @ Table[ND[invFreeGreenFunc[k], k, poles[[i]]], {i, Length @ poles}]);
		TensorMultiply[-I/2 HintqHintq\[LetterSpace]kpole\[Nu]\[Nu]p00, poleWeights, {{1, 1}}]
		}
	]
	
];


(* ::Section:: *)
(*7. Hubbard U calculation*)


(* ::Text:: *)
(*	Main functions of this package used to calculate the Hubbard on-site interaction energy U*)


(* ::Text:: *)
(*Quasi-1D on-shell T matrix calculation*)


setupTMatrixOnShellQuasi1D[vup_, vdown_, nCMmax_, nRelmax_, nkReg_, nkPole_, pOnShell_, \[Omega]perp_, smax_]:=Module[
	{
		E0, \[CapitalLambda]N, \[CapitalPi]HighEnergyCorrection, kReg, dkReg, kRegM, dkRegM, nM, dnM, counterTerm, Hintq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2,
		\[Phi]eigs0up\[LetterSpace]k\[Nu], \[Phi]eigs0down\[LetterSpace]k\[Nu], freeAtomsGFuncqsgeq1\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2, HintqFreeAtomsGFuncqsgeq1Weighted\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2, 
		\[CapitalPi]NonPoleMatrixgeq1, freeAtomsGFuncq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2, HintqFreeAtomsGFuncqWeighted\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2, \[CapitalPi]NonPoleMatrixseq0,
		\[CapitalPi]PoleMatrix, \[CapitalPi]Matrix, \[Eta]eigs\[LetterSpace]\[Nu], \[Eta]ergs\[LetterSpace]\[Nu]n, Hintqp\[LetterSpace]\[Nu]00, TMatrixOnShellQuasi1D, lowEnergyGridSize, qtildeN,
		highEnergyIntegrals\[LetterSpace]N, \[Eta]ergsHighEnergyIntegrals\[LetterSpace]\[Nu]n, highEnergyIntegralsQuasi\[LetterSpace]N, \[Eta]ergsHighEnergyIntegralsQuasi\[LetterSpace]\[Nu]n,
		\[CapitalPi]HighEnergyCorrectionQuasi
	},
	
	lowEnergyGridSize = Ceiling[nCMmax/2] + nRelmax + 1;
	{kReg, dkReg} = setupKIntReg[nkReg];
	E0 = Ek0\[Sigma][vup, lowEnergyGridSize, - pOnShell] + Ek0\[Sigma][vdown, lowEnergyGridSize, pOnShell];
	
	
	{\[Eta]eigs\[LetterSpace]\[Nu], \[Eta]ergs\[LetterSpace]\[Nu]n} = Chop @ EigSystHc[vup, vdown, nCMmax, 0];
	\[Eta]ergs\[LetterSpace]\[Nu]n = SparseArray @ \[Eta]ergs\[LetterSpace]\[Nu]n;
	
	(*High energy correction s=0*)
	
	highEnergyIntegrals\[LetterSpace]N = Re@Table[
		\[CapitalLambda]N = If[EvenQ @ Nvar, \[Pi] (2nRelmax + 1), 2\[Pi] nRelmax];
		1/(\[Pi] Sqrt[4 E0 - (2\[Pi] Nvar)^2] + 10^-6) Log[-((Sqrt[4 E0 - (2\[Pi] Nvar)^2]- 2 \[CapitalLambda]N)/(Sqrt[4 E0 - (2\[Pi] Nvar)^2]+ 2 \[CapitalLambda]N)) + 10^-6], (*check, especially Re@*)
		{Nvar, -nCMmax, nCMmax}];
	\[Eta]ergsHighEnergyIntegrals\[LetterSpace]\[Nu]n = highEnergyIntegrals\[LetterSpace]N # &/@ \[Eta]ergs\[LetterSpace]\[Nu]n; (*check*)
	\[CapitalPi]HighEnergyCorrection = TensorMultiply[\[Eta]ergs\[LetterSpace]\[Nu]n, \[Eta]ergsHighEnergyIntegrals\[LetterSpace]\[Nu]n, {{2, 2}}];
	
	(*High energy correction s>=1*)
	
	highEnergyIntegralsQuasi\[LetterSpace]N = Table[
		qtildeN = 2\[Pi] Nvar;
		\[CapitalLambda]N = If[EvenQ @ Nvar, \[Pi] (2nRelmax + 1), 2\[Pi] nRelmax];
		1/\[Pi] NIntegrate[quasi1DSum[E0 - (2\[Pi] Nvar)^2/4 - k^2, -\[Omega]perp, 1, smax], {k, \[CapitalLambda]N, \[Infinity]}], (*check*)
		{Nvar, -nCMmax, nCMmax}];
		
	\[Eta]ergsHighEnergyIntegralsQuasi\[LetterSpace]\[Nu]n = highEnergyIntegralsQuasi\[LetterSpace]N # &/@ \[Eta]ergs\[LetterSpace]\[Nu]n; (*check*)
		
	\[CapitalPi]HighEnergyCorrectionQuasi = TensorMultiply[\[Eta]ergs\[LetterSpace]\[Nu]n, \[Eta]ergsHighEnergyIntegralsQuasi\[LetterSpace]\[Nu]n, {{2, 2}}];
	
	counterTerm = 1/(2\[Pi]) NIntegrate[quasi1DSum[k^2, \[Omega]perp, 1, smax], {k, -\[Infinity], \[Infinity]}];
	
	(*s\[GreaterEqual]1 term*)
	
	{Hintq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2, \[Phi]eigs0up\[LetterSpace]k\[Nu], \[Phi]eigs0down\[LetterSpace]k\[Nu]} = buildHintTensor[vup, vdown, kReg, nCMmax, nRelmax];	
	freeAtomsGFuncqsgeq1\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2 = Table[
		ConstantArray[Outer[quasi1DSum[E0 - #1 - #2, -\[Omega]perp, 1, smax]&, \[Phi]eigs0up\[LetterSpace]k\[Nu][[ki]], \[Phi]eigs0down\[LetterSpace]k\[Nu][[ki]]], 2nCMmax + 1],
		{ki, Length @ kReg}];
	HintqFreeAtomsGFuncqsgeq1Weighted\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2 = dkReg/(2\[Pi]) * Hintq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2 * freeAtomsGFuncqsgeq1\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2; 
	
	\[CapitalPi]NonPoleMatrixgeq1 = TensorMultiply[
		Hintq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2,
		HintqFreeAtomsGFuncqsgeq1Weighted\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2,
		{{1, 1}, {3, 3}, {4, 4}}
	] + counterTerm * IdentityMatrix[2nCMmax + 1] + \[CapitalPi]HighEnergyCorrectionQuasi;
	
	
	(*s=0 term*)
	
	Hintq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2[[;;, ;;, 1, 1]] = 0;(*remove pole contribution*)
	freeAtomsGFuncq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2 = Table[
		ConstantArray[Outer[1 / (E0 - #1 - #2)&, \[Phi]eigs0up\[LetterSpace]k\[Nu][[ki]], \[Phi]eigs0down\[LetterSpace]k\[Nu][[ki]]], 2nCMmax + 1],
		{ki, Length @ kReg}];
	HintqFreeAtomsGFuncqWeighted\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2 = dkReg/(2\[Pi]) * Hintq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2 * freeAtomsGFuncq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2; 
	
	\[CapitalPi]NonPoleMatrixseq0 = TensorMultiply[
		Hintq\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2,
		HintqFreeAtomsGFuncqWeighted\[LetterSpace]k\[Nu]\[Nu]1\[Nu]2,
		{{1, 1}, {3, 3}, {4, 4}}
	] + \[CapitalPi]HighEnergyCorrection (** IdentityMatrix[2nCMmax+1]*); 
	
	
	(*POLE CONTRIBUTION*)
	\[CapitalPi]PoleMatrix = \[CapitalPi]PoleComponent[vup, vdown, pOnShell, nkPole, nCMmax, nRelmax];
	
	(*T Matrix*)
	\[CapitalPi]Matrix = \[CapitalPi]NonPoleMatrixseq0 + \[CapitalPi]NonPoleMatrixgeq1 + \[CapitalPi]PoleMatrix;
		
	Hintqp\[LetterSpace]\[Nu]00 = buildHintTensor\[Nu]1\[Nu]2eq00[vup, vdown, pOnShell, nCMmax, nRelmax, \[Eta]ergs\[LetterSpace]\[Nu]n];


	TMatrixOnShellQuasi1D[ainv_] := Hintqp\[LetterSpace]\[Nu]00 . Inverse[1/(-2 ainv) * IdentityMatrix[Length @ \[CapitalPi]Matrix] - \[CapitalPi]Matrix] . Hintqp\[LetterSpace]\[Nu]00;
	TMatrixOnShellQuasi1D
	
]


(* ::Text:: *)
(*Quasi-1D Hubbard U calculation*)


setupHubbardUQuasi1DAllParams[vup_, vdown_, nCMmax_, nRelmax_, nkReg_, nkPole_, pOnShell_, \[Omega]perp_, smax_] := Module[
	{
		TMatrixOnShellQuasi1D, tup, tdown, hubbardTMatrixInt, effMassHubbard, hubbardUQuasi1D,
		lowEnergyGridSize
	},
	
	TMatrixOnShellQuasi1D = setupTMatrixOnShellQuasi1D[vup, vdown, nCMmax, nRelmax, nkReg, nkPole, pOnShell, \[Omega]perp, smax];
	
	lowEnergyGridSize = Ceiling[nCMmax/2] + nRelmax + 1;
	tup = -1 / (2\[Pi]) NIntegrate[Ek0\[Sigma][vup, lowEnergyGridSize, k] Cos[k], {k, -\[Pi], \[Pi]}];
	tdown = -1 / (2\[Pi]) NIntegrate[Ek0\[Sigma][vdown, lowEnergyGridSize, k] Cos[k], {k, -\[Pi], \[Pi]}];
	
	hubbardTMatrixInt = HubbardTMatrixIntegral[tup, tdown, pOnShell];
	
	effMassHubbard = Im[hubbardTMatrixInt] / Im[TMatrixOnShellQuasi1D[-20]^-1];
	
	hubbardUQuasi1D[ainv_] := (effMassHubbard Re[TMatrixOnShellQuasi1D[ainv]^-1])^-1; (*rename to 1D*)
	
	hubbardUQuasi1D
]


(* ::Text:: *)
(*Quasi-1D Hubbard U calculation with automated convergence*)


setupHubbardUQuasi1D[vup_, vdown_, pOnShell_, \[Omega]perp_, convParam_] := setupHubbardUQuasi1DAllParams[vup, vdown, convParam * 5 + 5, convParam * 5 + 12, convParam * 3 + 8, 5 * convParam + 20, pOnShell, \[Omega]perp, 100 * convParam + 200]
