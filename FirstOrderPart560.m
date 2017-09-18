(* :Title: Cosmology36/NumericalPerturbation/FirstOrderPart560.m *)
(* :Notebook: Cosmology 36, p.31 *)

(* :Author: Yu-Hsiang Lin *)
(* :Date: Thu, Apr 16, 2015 *)

(* :Script Version: 55.0 *)
(* :Mathematica Version: 9.0 *)

(* :Context: Global` *)

(*
:Requirements:

Cosmology30/NumericalPerturbation/solveHybridPerturbation203.m
*)

(*
:Summary:

Find the spectrum of the curvature perturbations.
*)

(*
:Keywords:


*)

(*
:Discussion:

Version 32.0:
	Automate.
	This is the first version of the FirstOrderPart-series.

Version 33.0:
	Parallelize.

Version 34.0:
	Parameter search.
	For each lambdaPrime, fit A_s and n_s by m and initialPhiBar.

Version 34.1:
	Check lambdaPrime = 0 and lambda = 0, see how it goes.

Version 34.2:
	Debug.
	Keep two-field potential, but set psi to 0 everywhere.

Version 34.3:
	Totally single-field.

Version 35.0:
	Fix the phase bug.

Version 36.0:
	Discard interpolation for the coefficients.
	=> Idiot. This of course does not work.
	First, it has nothing to do with k-dependent phenomena.
	Second, it eats up all the memory still.

Version 37.0:
	I found that changing initial time seems to be helpful.
	Let me see unequal-time spectrum.

Version 38.0:
	Now the bug is confirmed to be numerical issue.
	We should reduce the time the mode being numerically evaluated within sub-horizon regime.
	Now go back to two-field.

Version 39.0:
	Push the numerical evaluation starting time to the horizon crossing.
	Increase the resolution in the interesting k range.

Version 39.1:
	Find out who is responsible for the peak.

Version 40.0:
	Problem identified: Don't use kSeparate!

Version 40.1:
	Find parameters.

Version 40.2:
	Find parameters.

Version 40.3:
	Find parameters.

Version 40.4:
	Find long-tail oscillations.

Version 40.5:
	Forget about long-tail. Keep searching parameters.

Version 40.6:
	Keep searching.

Version 40.7:
	Searching.

Version 41.0:
	Change the potential into "m^2 phi^2 / 2 + M^2 psi^2 / 2".

Version 42.0:
	Searching. Fine-tune.

Version 42.1:
	Searching. Fine-tune.

Version 43.0:
	Add e-folds after symmetry breaking.

Version 44.0:
	Try again "m^2 phi^2 / 2 + M^2 psi^2 / 2".

Version 45.0:
	Try \psi^4.

Version 46.0:
	Follow version 42.1, check calculation details.

Version 47.0:
	Kinetic domination. Single filed.

Version 47.1:
	All start from N = 0.

Version 47.2:
	No Hankel function. All with exponential initial condition.

Version 48.0:
	Going back to Version 40.7.
	Now my goal is to test for the stability.
	We can change the number of k modes sampled, the nBegin/nEnd.
	In this version I first change the nBegin and nEnd with lower k resolution.
	Conclusion: N=2 shift is probably a stable point.

Version 48.1:
	With N=2, change number of k modes sampled.
	Conclusion: 2x sampling number of k would be adequate to tell whether there is a numerical discontinuity.

Version 48.2:
	Try 2x sampling number + N=5 shift. See whether the discontinuity disappear.
	Conclusion: Yes, this prescription is consistent with the previous results and removes the numerical discontinuity.

Version 49.0:
	Only improve the middle part of k range sampling number by 2x. Use N=5 shift.
	Scan lambdaPrime as in Version 29.X.

Version 50.0:
	Scan initialPhiBar as in Version 28.X.

Version 51.0:
	Scan initialPsiBar as in Version 27.X.

Version 52.0:
	Use (coupling term + m^2 phi^2 term), not hybrid potential.
	Scan lambdaPrime.

Version 53.0:
	Use (coupling term + m^2 phi^2 term), not hybrid potential.
	Scan initialPhiBar.

Version 54.0:
	Use (coupling term + m^2 phi^2 term), not hybrid potential.
	Scan initialPsiBar.

Version 55.0:
	Use (coupling term + m^2 phi^2 term), not hybrid potential.
	Scan initialPhiBar. Just for initialPhiBar = 3.63.

Version 55.0:
	Use (coupling term + m^2 phi^2 term), not hybrid potential.
	Scan initialPsiBar. Just for initialPsiBar = 1.64.
*)



(* //////////////////////////////////////////////
(* The package loading has been moved to the Main code. *)
(* ------------------------------------------- *)
(* Load packages *)
(* ------------------------------------------- *)

Needs["solveHybridPerturbation201`"];
////////////////////////////////////////////// *)


(* ------------------------------------------- *)
(* Specify the range of comoving wavenumber k's we are going to solve *)
(* ------------------------------------------- *)

(* The kinetic energy, Hubble parameter, and Hubble radius *)

XFound[lna_] = 
	1/2 *
	(
		8*Pi * potential[phiBarFound[lna], psiBarFound[lna]] *
		( phiBarFound'[lna]^2 + psiBarFound'[lna]^2 ) /
		( 3 - 4*Pi * ( phiBarFound'[lna]^2 + psiBarFound'[lna]^2 ) )
	);

H[lna_] =
	(8*Pi/3)^(1/2) *
	( XFound[lna] + potential[phiBarFound[lna], psiBarFound[lna]] )^(1/2);

figrHLog =
	Plot[
		Log[1/H[n]],
		{n, 0, n34},
		
		PlotStyle -> RGBColor[1, 0, 0],
		PlotRange -> All,
		Axes -> False, Frame -> True,
		FrameLabel -> {"N", "ln(1/H)"},
		ImageSize -> 350,
		PlotPoints -> 150,
		AspectRatio -> Automatic
	];

Print[
	"\n",
	"figrHLog = \n",
	figrHLog
];



(* ------------------------------------------- *)

(* Expressing the comoving wavenumber k in terms of the Hubble radius today *)

(* The number of e-folds from the end of reheating to the present day *)

LCDMEFold = 65.8959;



(* The number of e-folds from the symmetry breaking to the end of reheating *)

endingEFold = 0;



(* Added by the e-folds from t_i to the symmetry breaking, n34, we obtained the total number of e-folds from t_i to the present day *)

initialA = E^( -( n34 + endingEFold + LCDMEFold ) );



(* The Hubble constant today (km / s / Mpc), by Planck 2013 *)

HubbleParameterToday = 67.11;



(* The comoving wavenumber k is how many times of the Hubble radius today *)

kInRH[n_] = 2*Pi/(n*c*Mpc/(HubbleParameterToday*1000*PlanckLength));



(* ------------------------------------------- *)

(* Find kSeparate before which the modes exit the pre-inflation horizon
and after which the modes exit the inflation horizon *)

lnrHPrime[n_] = D[ Log[ 1 / H[n] ], n ];



(* There's more than one possible solutions... *)

nTrialLow = 15;
nTrialHigh = 20;

nSeparateFoundList =
	Table[
		n /. FindRoot[
			lnrHPrime[n] == 1,
			
			{n, nTrial}
		],
		
		{nTrial, nTrialLow, nTrialHigh}
	];



(* ... except the nonsense ones... *)

nSeparateFoundListFiltered = {};

Do[
	If[ ( nSeparateFoundList[[index]] > nTrialLow )  &&  ( nSeparateFoundList[[index]] < nTrialHigh ),
		nSeparateFoundListFiltered =
			Join[ nSeparateFoundListFiltered, {nSeparateFoundList[[index]]} ]
	],
	
	{index, 1, Length[nSeparateFoundList]}
];



(* ... we want the smallest *)

nSeparateFound = Union[nSeparateFoundListFiltered][[1]];

Print[
	"\n",
	"nSeparateFound = ", nSeparateFound
];

kSeparate =
	initialA * Exp[nSeparateFound] * H[nSeparateFound];



(* ------------------------------------------- *)

(* The Planck pivot k---"k_0" *)

kPivot = 0.05 * PlanckLength / Mpc;



(* The range of k we are going to solve *)

kListSmall =
	Table[
		10^logn * PlanckLength / Mpc,
		
		{
			logn,
			-5.6,
			Log[10, kSeparate * Mpc / PlanckLength] - 1 - 0.025,
			0.05
		}   (* by experience *)
	];

kListTransition =
	Table[
		10^logn * PlanckLength / Mpc,
		
		{
			logn,
			Log[10, kSeparate * Mpc / PlanckLength] - 1 - 0.01,
			Log[10, kSeparate * Mpc / PlanckLength] + 1,
			0.0125
			(* original interval 0.025 *)
		}   (* by experience *)
	];

kListLarge =
	Table[
		10^logn * PlanckLength / Mpc,
		
		{
			logn,
			Log[10, kSeparate * Mpc / PlanckLength] + 1 + 0.025,
			0.6,
			0.1
		}   (* by experience *)
	];



kList =
	Join[
		kListSmall, kListTransition, kListLarge
	];



Print[
	"\n",
	"The length of kListSmall = ", Length[kListSmall], "\n",
	"The length of kListTransition = ", Length[kListTransition], "\n",
	"The length of kListLarge = ", Length[kListLarge], "\n",
	"The total length of kList's = ", Length[kList]
];



(* ------------------------------------------- *)

(* { k1, k2, ..., kIndexSmallEnd} => small k, long wavelength, exit horizon earlier, use nBeginSmall
   { kIndexSmallEnd + 1, ..., Length[kList]} => large k, short wavelength, exit horizon later, use nBeginLarge *)

kIndexTransitionSmallEnd = Length[kListTransition];

While[ ( kListTransition[[kIndexTransitionSmallEnd]] > kSeparate ) && ( kIndexTransitionSmallEnd >= 2 ),
	kIndexTransitionSmallEnd = kIndexTransitionSmallEnd - 1;
];



(* Set the N from which I begin the numerical solution *)

nBeginListSmall =
	Table[
		( n - 4 ) /. FindRoot[
			Log[kListSmall[[kIndex]]] == Log[initialA] + n + Log[H[n]],
			{n, nSeparateFound - 1-1}   (* From experience *)
		],
		
		{kIndex, 1, Length[kListSmall]}
	];



nBeginListTransitionSmall =
	Table[
		( n - 4 ) /. FindRoot[
			Log[kListTransition[[kIndex]]] == Log[initialA] + n + Log[H[n]],
			{n, nSeparateFound - 1-1}   (* From experience *)
		],
		
		{kIndex, 1, kIndexTransitionSmallEnd}
	];

nBeginListTransitionLarge =
	Table[
		( n - 4 ) /. FindRoot[
			Log[kListTransition[[kIndex]]] == Log[initialA] + n + Log[H[n]],
			{n, nSeparateFound + 3-1}   (* From experience *)
		],
		
		{kIndex, kIndexTransitionSmallEnd + 1, Length[kListTransition]}
	];

nBeginListTransition =
	Join[
		nBeginListTransitionSmall,
		nBeginListTransitionLarge
	];



nBeginListLarge =
	Table[
		( n - 4 ) /. FindRoot[
			Log[kListLarge[[kIndex]]] == Log[initialA] + n + Log[H[n]],
			{n, nSeparateFound + 3-1}   (* From experience *)
		],
		
		{kIndex, 1, Length[kListLarge]}
	];



nBeginList =
	Join[
		nBeginListSmall,
		nBeginListTransition,
		nBeginListLarge
	];



(* Here I want to check that all nBegin are larger than zero *)

Print[
	"\n",
	"nBeginList = ", nBeginList
];



(* ------------------------------------------- *)

(* Plot the evolution of the physical wavelength against the number of e-fold *)

lnLambdaPhys[k_, n_] =
	Log[initialA / k] + n;

nbegin = 0;
nend = n34;



figModeListSmall =
	Table[
		Plot[
			lnLambdaPhys[kListSmall[[kIndex]], n],
			{n, nbegin, nend},
			
			PlotRange -> All,
			Axes -> False, Frame -> True, FrameTicks -> All,
			PlotStyle -> {ColorData["Rainbow"][0.55]}
		],
		
		{kIndex, 1, Length[kListSmall]}
	];

figNBeginListSmall =
	Table[
		Graphics[{
			PointSize[0.005], ColorData["Rainbow"][0.55],
			Point[{
				nBeginListSmall[[kIndex]],
				lnLambdaPhys[kListSmall[[kIndex]], nBeginListSmall[[kIndex]]]
			}]
		}],
		
		{kIndex, 1, Length[kListSmall]}
	];



figModeListTransitionSmall =
	Table[
		Plot[
			lnLambdaPhys[kListTransition[[kIndex]], n],
			{n, nbegin, nend},
			
			PlotRange -> All,
			Axes -> False, Frame -> True, FrameTicks -> All,
			PlotStyle -> {ColorData["Rainbow"][0.55]}
		],
		
		{kIndex, 1, kIndexTransitionSmallEnd}
	];

figNBeginListTransitionSmall =
	Table[
		Graphics[{
			PointSize[0.005], ColorData["Rainbow"][0.55],
			Point[{
				nBeginListTransition[[kIndex]],
				lnLambdaPhys[kListTransition[[kIndex]], nBeginListTransition[[kIndex]]]
			}]
		}],
		
		{kIndex, 1, kIndexTransitionSmallEnd}
	];

figModeListTransitionLarge =
	Table[
		Plot[
			lnLambdaPhys[kListTransition[[kIndex]], n],
			{n, nbegin, nend},
			
			PlotRange -> All,
			Axes -> False, Frame -> True, FrameTicks -> All,
			PlotStyle -> {RGBColor[0, 0, 0]}
		],
		
		{kIndex, kIndexTransitionSmallEnd + 1, Length[kListTransition]}
	];

figNBeginListTransitionLarge =
	Table[
		Graphics[{
			PointSize[0.005], RGBColor[0, 0, 0],
			Point[{
				nBeginListTransition[[kIndex]],
				lnLambdaPhys[kListTransition[[kIndex]], nBeginListTransition[[kIndex]]]
			}]
		}],
		
		{kIndex, kIndexTransitionSmallEnd + 1, Length[kListTransition]}
	];



figModeListLarge =
	Table[
		Plot[
			lnLambdaPhys[kListLarge[[kIndex]], n],
			{n, nbegin, nend},
			
			PlotRange -> All,
			Axes -> False, Frame -> True, FrameTicks -> All,
			PlotStyle -> {RGBColor[0, 0, 0]}
		],
		
		{kIndex, 1, Length[kListLarge]}
	];

figNBeginListLarge =
	Table[
		Graphics[{
			PointSize[0.005], RGBColor[0, 0, 0],
			Point[{
				nBeginListLarge[[kIndex]],
				lnLambdaPhys[kListLarge[[kIndex]], nBeginListLarge[[kIndex]]]
			}]
		}],
		
		{kIndex, 1, Length[kListLarge]}
	];



figHorizonExit = 
	Show[
		figrHLog,
		figNBeginListSmall,
		figNBeginListTransitionSmall,
		figNBeginListTransitionLarge,
		figNBeginListLarge,
		
		FrameLabel -> {"N", "ln(1/H)"},
		ImageSize -> 1000, PlotRange -> {{0, 35}, {-1, 15}},
		AspectRatio -> 1 / GoldenRatio
	];


(*
figHorizonExit = 
	Show[
		figrHLog,
		figModeListSmall, figNBeginListSmall,
		figModeListTransitionSmall, figNBeginListTransitionSmall,
		figModeListTransitionLarge, figNBeginListTransitionLarge,
		figModeListLarge, figNBeginListLarge,
		
		FrameLabel -> {"N", "ln(1/H)"},
		ImageSize -> 1000, PlotRange -> {{0, 35}, {-1, 15}}
	];
*)

Print[
	"\n",
	"figHorizonExit = \n",
	figHorizonExit
];



SetDirectory[outputFolder];

Export[
	"figHorizonExit" <>
		IntegerString[lambdaPrimeIndex, 10, 2] <>
		IntegerString[mIndex, 10, 2] <>
		IntegerString[initialPsiBarIndex, 10, 2] <>
		".pdf",
	figHorizonExit
];

ResetDirectory[];



(* ------------------------------------------- *)
(* Numerical solution *)
(* ------------------------------------------- *)

(* Calculate the conformal time *)

solEta =
	NDSolve[
		{
			eta'[lna] == 1 / ( initialA * Exp[lna] * H[lna] ),
			
			eta[0] == 0
		},
		
		eta,
		{lna, 0, n34}
	];

etaFound[lna_] = Re[eta[lna]] /. solEta[[1]];   (* It happens that there is a null imaginary part in the output of NDSolve[] *)



(* ------------------------------------------- *)

(* The initial phase *)

(* In our vanilla model, it is uniformly zero *)

deltaPhiInitialPhaseList =
	Table[
		0,
		
		{kIndex, 1, Length[kList]}
	];

deltaPsiInitialPhaseList =
	Table[
		0,
		
		{kIndex, 1, Length[kList]}
	];



(* ------------------------------------------- *)

(* From N = 0 to nBegin, we have good approximation for sub-horizon modes *)

deltaPhiBeginList =
	Table[
		1 / ( (2*Pi)^(3/2) * ( 2 * kList[[kIndex]] )^(1/2) ) *
		1 / ( initialA * Exp[  nBeginList[[kIndex]]  ]) *
		Exp[ -I * kList[[kIndex]] * etaFound[nBeginList[[kIndex]]] ] *
		Exp[ I * deltaPhiInitialPhaseList[[kIndex]] ],
		
		{kIndex, 1, Length[kList]}
	];

deltaPhiPrimeBeginList =
	Table[
		-(
			1 + I * kList[[kIndex]] /
				(
					initialA * Exp[  nBeginList[[kIndex]]  ] *
					H[  nBeginList[[kIndex]]  ]
				)
		) * deltaPhiBeginList[[kIndex]],
		
		{kIndex, 1, Length[kList]}
	];

deltaPsiBeginList =
	Table[
		1 / ( (2*Pi)^(3/2) * ( 2 * kList[[kIndex]] )^(1/2) ) *
		1 / ( initialA * Exp[  nBeginList[[kIndex]]  ]) *
		Exp[ -I * kList[[kIndex]] * etaFound[nBeginList[[kIndex]]] ] *
		Exp[ I * deltaPsiInitialPhaseList[[kIndex]] ],
		
		{kIndex, 1, Length[kList]}
	];

deltaPsiPrimeBeginList =
	Table[
		-(
			1 + I * kList[[kIndex]] /
				(
					initialA * Exp[  nBeginList[[kIndex]]  ] *
					H[  nBeginList[[kIndex]]  ]
				)
		) * deltaPsiBeginList[[kIndex]],
		
		{kIndex, 1, Length[kList]}
	];

PsiBeginList =
	Table[
		1 /
		(
			2 * XFound[  nBeginList[[kIndex]]  ] - 
			1 / (4*Pi) * kList[[kIndex]]^2 / initialA^2 * Exp[ - 2 * nBeginList[[kIndex]] ]
		) *
		(
			H[  nBeginList[[kIndex]]  ]^2 * 
				(
					phiBarFound'[nBeginList[[kIndex]]] * deltaPhiPrimeBeginList[[kIndex]] +
					psiBarFound'[nBeginList[[kIndex]]] * deltaPsiPrimeBeginList[[kIndex]]
				) +
			(
				3 * H[  nBeginList[[kIndex]]  ]^2 * phiBarFound'[nBeginList[[kIndex]]] +
				(  D[potential[phiBar, psiBar], phiBar] /.
					{phiBar -> phiBarFound[nBeginList[[kIndex]]], psiBar -> psiBarFound[nBeginList[[kIndex]]]}  )
			) * deltaPhiBeginList[[kIndex]] +
			(
				3 * H[  nBeginList[[kIndex]]  ]^2 * psiBarFound'[nBeginList[[kIndex]]] +
				(  D[potential[phiBar, psiBar], psiBar] /.
					{phiBar -> phiBarFound[nBeginList[[kIndex]]], psiBar -> psiBarFound[nBeginList[[kIndex]]]}  )
			) * deltaPsiBeginList[[kIndex]]
		),
		
		{kIndex, 1, Length[kList]}
	];



(* ------------------------------------------- *)

(* Set up parallel calculation *)

(* Launch kernels *)

kernelList = LaunchKernels[];

Print[
	"\n",
	"We have kernels: \n",
	kernelList
];



(* The value we are going to record the steady value *)

nSteady = 40;

nSolveEnd = nSteady + 0.1;



(* Distribute variables into parallel kernels *)

DistributeDefinitions[
	solvePerturbation,
	potential,
	phiBarFound, psiBarFound,
	deltaPhiBeginList, deltaPhiPrimeBeginList, 
	deltaPsiBeginList, deltaPsiPrimeBeginList, 
	PsiBeginList,
	initialA, kList, nBeginList, nSolveEnd
];



(* Make space to store results in the kernels *)

ParallelEvaluate[
	solPerturbationList =
		Table[
			0,
			
			{kIndex, 1, Length[kList]}
		];
];



(* ------------------------------------------- *)

(* From N = nBegin to n34, we solve numerically *)

(* Parallel calculation *)

timeBegin = DateString[];

Print[
	"\n",
	"Parallel solve begins at ", timeBegin
];



ParallelDo[
	solPerturbationList[[kIndex]] =
		solvePerturbation[
			potential, phiBarFound, psiBarFound,
			deltaPhiBeginList[[kIndex]], deltaPhiPrimeBeginList[[kIndex]], 
			deltaPsiBeginList[[kIndex]], deltaPsiPrimeBeginList[[kIndex]], 
			PsiBeginList[[kIndex]],
			initialA, kList[[kIndex]], nBeginList[[kIndex]], nSolveEnd
		];,
	
	{kIndex, 1, Length[kList]}
];



timeEnd = DateString[];

Print[
	"\n",
	"Parallel solve ends at ", timeEnd
];

timeUsed =
	DateDifference[
		timeBegin, timeEnd,
		"Minute"
	][[1]];

Print[
	"\n",
	"It takes ", timeUsed, " minutes"
];



(* ------------------------------------------- *)

(* Collect the parallel results *)

solPerturbationList =
	Apply[Plus, ParallelEvaluate[solPerturbationList]];
(* //////////////////////////////////////////////
Print[
	"\n",
	"solPerturbationList = \n",
	solPerturbationList
];
////////////////////////////////////////////// *)


(* ------------------------------------------- *)

(* Close the kernels *)

CloseKernels[];



(* ------------------------------------------- *)
(* The spectrum of the curvature perturbations *)
(* ------------------------------------------- *)

(* Record value at N = 60, arbitrarily chosen steady value *)

(* nSteady = 30; *)

chiSteadyList =
	Table[
		( solPerturbationList[[kIndex, 5]][nSteady] + I * solPerturbationList[[kIndex, 6]][nSteady] ) +
		4*Pi/3 *
			( 1 + potential[phiBarFound[nSteady], psiBarFound[nSteady]] / XFound[nSteady] ) *
			(
				phiBarFound'[nSteady] *
					( solPerturbationList[[kIndex, 1]][nSteady] + I * solPerturbationList[[kIndex, 2]][nSteady] ) +
				psiBarFound'[nSteady] *
					( solPerturbationList[[kIndex, 3]][nSteady] + I * solPerturbationList[[kIndex, 4]][nSteady] )
			),
		
		{kIndex, 1, Length[kList]}
	];



(* ------------------------------------------- *)

(* Clean up solPerturbationList in the main kernel. *)

Clear[solPerturbationList];

Remove["Global`solPerturbationList"];



(* ------------------------------------------- *)

(* The spectrum of the curvature perturbations *)

chiSpectrumList =
	Table[
		{
			Log[10, kList[[kIndex]] * Mpc / PlanckLength],
			4*Pi * kList[[kIndex]]^3 * Abs[chiSteadyList[[kIndex]]]^2
		},
		
		{kIndex, 1, Length[kList]}
	];

chiSpectrumFunction = 
	Interpolation[chiSpectrumList, InterpolationOrder -> 1];



(* ------------------------------------------- *)

(* Calculate the amplitude and the spectral index at the pivot k *)

as =
	chiSpectrumFunction[  Log[10, kPivot * Mpc / PlanckLength]  ];

ns =
	1 +
	D[Log[10, chiSpectrumFunction[logkVar]], logkVar] /.
		{logkVar -> Log[10, kPivot * Mpc / PlanckLength]};

Print[
	"\n",
	"kPivot = ", kPivot * Mpc / PlanckLength, " \!\(\*SuperscriptBox[\(Mpc\), \(-1\)]\)", "\n",
	"\!\(\*SubscriptBox[\(A\), \(s\)]\) = ", as, "\n",
	"\!\(\*SubscriptBox[\(n\), \(s\)]\) = ", ns
];



(* ------------------------------------------- *)

data =
	Table[
		{
			Log[10, kList[[kIndex]] * Mpc / PlanckLength],
			Log[10, chiSpectrumList[[kIndex, 2]]]
		},
		
		{kIndex, 1, Length[kList]}
	];



(* ------------------------------------------- *)

(* Plot the spectrum of the curvature perturbations *)

figChiSpectrum =
	ListPlot[
		data,
		
		Axes -> False, Frame -> True, FrameTicks -> All,
		FrameLabel ->
			{
				"\!\(\*SubscriptBox[\(log\), \(10\)]\)(k/\!\(\*SuperscriptBox[\(Mpc\), \(-1\)]\))",
				"\!\(\*SubscriptBox[\(log\), \(10\)]\)(4\!\(\*SuperscriptBox[\(\[Pi]k\), \(3\)]\)|\[Chi]\!\(\*SuperscriptBox[\(|\), \(2\)]\))"
			},
		ImageSize -> 350, PlotRange -> All, Joined -> True,
		InterpolationOrder -> 1, Mesh -> Full
	];

figChiSpectrumZoomedIn = 
	ListPlot[
		Take[data, -Round[  0.6 * Length[kList]  ]],
		
		Axes -> False, Frame -> True, 
		FrameTicks -> All, 
		FrameLabel ->
		{
			"\!\(\*SubscriptBox[\(log\), \(10\)]\)(k/\!\(\*SuperscriptBox[\(Mpc\), \(-1\)]\))",
			"\!\(\*SubscriptBox[\(log\), \(10\)]\)(4\!\(\*SuperscriptBox[\(\[Pi]k\), \(3\)]\)|\[Chi]\!\(\*SuperscriptBox[\(|\), \(2\)]\))"
		},
		ImageSize -> 700, PlotRange -> All, Joined -> True, 
		InterpolationOrder -> 1, Mesh -> Full
	];

Print[
	"\n",
	"figChiSpectrum = \n",
	figChiSpectrum
];

Print[
	"\n",
	"figChiSpectrumZoomedIn = \n",
	figChiSpectrumZoomedIn
];



(* ------------------------------------------- *)

(* Export chiSpectrumList and the figures *)

SetDirectory[outputFolder];

DumpSave[
	"chiSpectrumList" <>
		IntegerString[lambdaPrimeIndex, 10, 2] <>
		IntegerString[mIndex, 10, 2] <>
		IntegerString[initialPsiBarIndex, 10, 2] <>
		".mx",
	chiSpectrumList
];

Export[
	"figChiSpectrum" <>
		IntegerString[lambdaPrimeIndex, 10, 2] <>
		IntegerString[mIndex, 10, 2] <>
		IntegerString[initialPsiBarIndex, 10, 2] <>
		".pdf",
	
	figChiSpectrum
];

Export[
	"figChiSpectrumZoomedIn" <>
		IntegerString[lambdaPrimeIndex, 10, 2] <>
		IntegerString[mIndex, 10, 2] <>
		IntegerString[initialPsiBarIndex, 10, 2] <>
		".pdf",
	
	figChiSpectrumZoomedIn
];

ResetDirectory[];



(* ------------------------------------------- *)
(* Make the CAMB input files *)
(* ------------------------------------------- *)

{xList, yList} = Transpose[data];



(* xList *)

Do[
	If[ xList[[index]] < 0,
		(
			xList[[index]] =
				StringJoin[
					Insert[
						Insert[
							ToString /@ RealDigits[  xList[[index]], 10, 8, 0  ][[1]],
							".", 2
						],
						"-", 1
					]
				];
		),
		
		(* xList[[index]] >= 0 *)
		(
			xList[[index]] =
				StringJoin[
					Insert[
						ToString /@ RealDigits[  xList[[index]], 10, 8, 0  ][[1]],
						".", 2
					]
				];
		)
	],
	
	{index, 1, Length[xList]}
];

posList = Table[{k}, {k, 2, Length[xList] + 1}];

xList = Insert[xList, "\n", posList];

xList = StringJoin[xList];



(* yList *)

yList = 10^yList;

Do[
	yList[[index]] =
		StringJoin[
			StringJoin[
				Insert[
					ToString /@ RealDigits[  yList[[index]], 10, 8, Floor[Log[10, yList[[index]]]]  ][[1]],
					".",2
				]
			],
			"E",
			ToString[  Floor[Log[10, yList[[index]]]]  ]
		];,
	
	{index, 1, Length[yList]}
];

posList = Table[{k}, {k, 2, Length[yList] + 1}];

yList = Insert[yList, "\n", posList];

yList = StringJoin[yList];



SetDirectory[outputFolder];

Export[
	"logkListFortran" <>
		IntegerString[lambdaPrimeIndex, 10, 2] <>
		IntegerString[mIndex, 10, 2] <>
		IntegerString[initialPsiBarIndex, 10, 2] <>
		".txt",
	xList
];

Export[
	"powerListFortran" <>
		IntegerString[lambdaPrimeIndex, 10, 2] <>
		IntegerString[mIndex, 10, 2] <>
		IntegerString[initialPsiBarIndex, 10, 2] <>
		".txt",
	yList
];

ResetDirectory[];
