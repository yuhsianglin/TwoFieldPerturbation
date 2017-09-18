(* :Title: Cosmology30/Numerical perturbation/solveHybridPerturbation203.m *)
(* :Notebook: Cosmology 30, p.77 *)

(* :Author: Yu-Hsiang Lin *)
(* :Date: Sat, Jul 19, 2014 *)

(* :Package Version: 20.3 *)
(* :Mathematica Version: 9.0 *)

(* :Context: solveHybridPerturbation203` *)

(*
:Requirements:


*)

(*
:Summary:

Solve the evolution of perturbations in inflation era. The variable is not cosmic time, but ln(a/ai); i.e. the e-fold.

Using hybrid inflation model.

Note that all quantities in this program are in
	Planck units.
*)

(*
:Keywords:


*)

(*
:Discussion:

Version 1.1:
	Reduce MaxSteps. Find out what the reason that makes it eat up all the memory is.

Version 2.0:
	Supposedly it's because the using of r*E^(I*theta) formalism that makes the program too complicated.
	In this version, switch back to real-part/imaginary-part formalism.

Version 2.1:
	Pass the solution.

Version 3.0:
	Change the variable into conformal time, eta.
	Treat Psi (the metric perturbation in Newtonian gauge) as another function to be solved through the differential equation.

Version 4.0:
	Use N-formalism.
	Debugging.
	Found the bug that there are just too many InterpolatingFunction!
	We collect (by hand) them and make a single InterpolatingFunction for each coefficient, then everything is set.

Version 4.1:
	The clean-up version.

Version 4.2:
	Migrate to Mathematica 9.
	Minor changes.

Version 10.0:
	Correct some errors.
	Change the initial conditions.

Version 17.0:
	Solve for the horizon exit period.
	We can assign initial values and derivatives of deltaPhi and deltaPsi when calling solvePerturbation[] now.

Version 17.2:
	Add time record.

Version 18.0:
	Do calculation without Which[].

Version 19.0:
	Calculate Psi by differential equation (instead of constraint equation, like what was done).

Version 20.0:
	Version 19.0 works, so we are now moving on calculating the perturbations.

Version 20.1:
	Suppress intermediate output.

Version 20.2:
	Debug: Turn off the psi field.

Version 20.3:
	Plot phase of the modes as well.

Version 21.0:
	Discard interpolation for the coefficients.
*)



(* -------------------------------------------- *)
(* -------------------------------------------- *)

BeginPackage["solveHybridPerturbation203`"];



(* -------------------------------------------- *)

(* Functions to be exported *)

{
	solvePerturbation
};



(* -------------------------------------------- *)
(* -------------------------------------------- *)

Begin["`Private`"];



(* -------------------------------------------- *)

(*
Solve for perturbations deltaPhi and deltaPsi.

Psi (metric perturbation) is now determined by constraint equation.

Note that phiBar(lna) and psiBar(lna) are both in terms of ln(a/ai).
*)

solvePerturbation[potential_, phiBar_, psiBar_, initialDeltaPhi_, initialDeltaPhiPrime_, initialDeltaPsi_, initialDeltaPsiPrime_, initialPsi_, initialA_, k_, lnaBegin_, lnaEnd_]:=
Module[
	{
		lna,
		
		X,
		reDeltaPhi, imDeltaPhi,
		reDeltaPsi, imDeltaPsi,
		
		HSquare,
		
		(*b, bValueList, bFunction,*)
		
		differentialEquations,
		
		deltaPhiInitial,
		deltaPsiInitial,
		PsiInitial,
		
		reDeltaPhiInitial, imDeltaPhiInitial,
		reDeltaPsiInitial, imDeltaPsiInitial,
		rePsiInitial, imPsiInitial,
		deltaPhiPrimeInitial,
		deltaPsiPrimeInitial,
		reDeltaPhiPrimeInitial, imDeltaPhiPrimeInitial,
		reDeltaPsiPrimeInitial, imDeltaPsiPrimeInitial,
		initialConditions,
		
		sol,
		reDeltaPhiFound, imDeltaPhiFound,
		reDeltaPsiFound, imDeltaPsiFound,
		lnaEndGot,
		
		Mpc, PlanckLength,
		figDeltaPhiAmplitudeSquared,
		figDeltaPsiAmplitudeSquared,
		
		timeBegin, timeEnd,
		timePreparing, timeSolving
	},
	
	
	
	(* -------------------------------------------- *)
	
	(* Print solving message *)
	
	Mpc = 30.857*10^15*10^6;
	PlanckLength = ( (1.054572*10^-34) * (6.67428*10^-11) / (2.99792458*10^8)^3 )^(1/2);
	
	(* Print["\n", "Solving the mode with k = ", k * Mpc / PlanckLength, " \!\(\*SuperscriptBox[\(Mpc\), \(-1\)]\)"]; *)
	
	
	
	(* -------------------------------------------- *)
	
	(* Record time for preparing *)
	
	timeBegin = DateList[];
	
	
	
	(* -------------------------------------------- *)
	
	(* The kinetic energy and Hubble parameter squared *)
	
	X[lna_] =
		4 * Pi * potential[phiBar[lna], psiBar[lna]] * ( phiBar'[lna]^2 + psiBar'[lna]^2 ) /
			(  3 - 4 * Pi * ( phiBar'[lna]^2 + psiBar'[lna]^2 )  );
	
	HSquare[lna_] =
		(8*Pi/3) * ( X[lna] + potential[phiBar[lna], psiBar[lna]] );
	
	
	
	(* -------------------------------------------- *)
	
	(* Declare memory space *)
	
	Array[b, {3, 9}];
	
	Array[bValueList, {3, 9}];
	
	Array[bFunction, {3, 9}];
	
	
	
	(* -------------------------------------------- *)
	
	(*
		Coefficients of the differential equations:
		
		b11 deltaPhi'' + b12 deltaPhi' + b13 deltaPhi +
		b14 deltaPsi'' + b15 deltaPsi' + b16 deltaPsi +
		b17 Psi''      + b18 Psi'      + b19 Psi        = 0
	*)
	
	
	b[1, 1][lna_] = 1;
	
	b[1, 2][lna_] = 8 * Pi * potential[phiBar[lna], psiBar[lna]] / HSquare[lna];
	
	b[1, 3][lna_] = ( k^2 / initialA^2 * Exp[-2*lna] + D[ potential[phiBar[lna], psiBar[lna]], {phiBar[lna], 2} ] ) / HSquare[lna];
	
	b[1, 4][lna_] = 0;
	
	b[1, 5][lna_] = 0;
	
	b[1, 6][lna_] = ( D[ D[ potential[phiBarVar, psiBarVar], phiBarVar ], psiBarVar ] /. { phiBarVar -> phiBar[lna], psiBarVar -> psiBar[lna] } ) / HSquare[lna];
	
	b[1, 7][lna_] = 0;
	
	b[1, 8][lna_] = -4 * phiBar'[lna];
	
	b[1, 9][lna_] = 2 * ( D[ potential[phiBarVar, psiBarVar], phiBarVar ] /. { phiBarVar -> phiBar[lna], psiBarVar -> psiBar[lna] } ) / HSquare[lna];
	
	
	
	b[2, 1][lna_] = 0;
	
	b[2, 2][lna_] = 0;
	
	b[2, 3][lna_] = ( D[ D[ potential[phiBarVar, psiBarVar], phiBarVar ], psiBarVar ] /. { phiBarVar -> phiBar[lna], psiBarVar -> psiBar[lna] } ) / HSquare[lna];
	
	b[2, 4][lna_] = 1;
	
	b[2, 5][lna_] = 8 * Pi * potential[phiBar[lna], psiBar[lna]] / HSquare[lna];
	
	b[2, 6][lna_] = ( k^2 / initialA^2 * Exp[-2*lna] + ( D[ potential[phiBarVar, psiBarVar], {psiBarVar, 2} ] /. { phiBarVar -> phiBar[lna], psiBarVar -> psiBar[lna] } ) ) / HSquare[lna];
	
	b[2, 7][lna_] = 0;
	
	b[2, 8][lna_] = -4 * psiBar'[lna];
	
	b[2, 9][lna_] = 2 * ( D[ potential[phiBarVar, psiBarVar], psiBarVar ] /. { phiBarVar -> phiBar[lna], psiBarVar -> psiBar[lna] } ) / HSquare[lna];
	
	
	
	b[3, 1][lna_] = 0;
	
	b[3, 2][lna_] = 0;
	
	b[3, 3][lna_] = -4 * Pi * phiBar'[lna];
	
	b[3, 4][lna_] = 0;
	
	b[3, 5][lna_] = 0;
	
	b[3, 6][lna_] = -4 * Pi * psiBar'[lna];
	
	b[3, 7][lna_] = 0;
	
	b[3, 8][lna_] = 1;
	
	b[3, 9][lna_] = 1;
	
	
	
	(* -------------------------------------------- *)
	
	(* Make each coefficient a interpolating function *)
	
	Do[(
		bValueList[index1, index2] =
			Table[
				{lna, b[index1, index2][lna]},
				
				{lna, lnaBegin, lnaEnd, (lnaEnd - lnaBegin) / 400}
			];
		
		bFunction[index1, index2] = Interpolation[bValueList[index1, index2], InterpolationOrder -> 1];
		
		(*Print[
			Plot[bFunction[index1, index2][lna], {lna, lnaBegin, lnaEnd}, PlotRange -> All, Axes -> False, Frame -> True, FrameTicks -> All, ImageSize -> 350, PlotLabel -> "{" <> ToString[index1] <> ", " <> ToString[index2] <> "}"]
		];
	
		Print[
			Plot[Log[10, Abs[bFunction[index1, index2][lna]]], {lna, lnaBegin, lnaEnd}, PlotRange -> All, Axes -> False, Frame -> True, FrameTicks -> All, ImageSize -> 350, PlotLabel -> "{" <> ToString[index1] <> ", " <> ToString[index2] <> "}"]
		];*)
		),
		
		{index1, 1, 3},
		{index2, 1, 9}
	];
	
	
	
	(* -------------------------------------------- *)
	
	(* Now we are ready for writing down the differential equations *)
	
	(* Note that we turn off Psi''[lna] by hand to keep NDSolve[] away from viewing Psi[] as a second order differential equation *)
	
	differentialEquations =
	{
		bFunction[1, 1][lna] * reDeltaPhi''[lna] + bFunction[1, 2][lna] * reDeltaPhi'[lna] + bFunction[1, 3][lna] * reDeltaPhi[lna] + bFunction[1, 4][lna] * reDeltaPsi''[lna] + bFunction[1, 5][lna] * reDeltaPsi'[lna] + bFunction[1, 6][lna] * reDeltaPsi[lna](* + bFunction[1, 7][lna] * rePsi''[lna]*) + bFunction[1, 8][lna] * rePsi'[lna] + bFunction[1, 9][lna] * rePsi[lna]== 0,
		
		bFunction[2, 1][lna] * reDeltaPhi''[lna] + bFunction[2, 2][lna] * reDeltaPhi'[lna] + bFunction[2, 3][lna] * reDeltaPhi[lna] + bFunction[2, 4][lna] * reDeltaPsi''[lna] + bFunction[2, 5][lna] * reDeltaPsi'[lna] + bFunction[2, 6][lna] * reDeltaPsi[lna](* + bFunction[2, 7][lna] * rePsi''[lna]*) + bFunction[2, 8][lna] * rePsi'[lna] + bFunction[2, 9][lna] * rePsi[lna]== 0,
		
		bFunction[3, 1][lna] * reDeltaPhi''[lna] + bFunction[3, 2][lna] * reDeltaPhi'[lna] + bFunction[3, 3][lna] * reDeltaPhi[lna] + bFunction[3, 4][lna] * reDeltaPsi''[lna] + bFunction[3, 5][lna] * reDeltaPsi'[lna] + bFunction[3, 6][lna] * reDeltaPsi[lna](* + bFunction[3, 7][lna] * rePsi''[lna]*) + bFunction[3, 8][lna] * rePsi'[lna] + bFunction[3, 9][lna] * rePsi[lna]== 0,
		
		bFunction[1, 1][lna] * imDeltaPhi''[lna] + bFunction[1, 2][lna] * imDeltaPhi'[lna] + bFunction[1, 3][lna] * imDeltaPhi[lna] + bFunction[1, 4][lna] * imDeltaPsi''[lna] + bFunction[1, 5][lna] * imDeltaPsi'[lna] + bFunction[1, 6][lna] * imDeltaPsi[lna](* + bFunction[1, 7][lna] * imPsi''[lna]*) + bFunction[1, 8][lna] * imPsi'[lna] + bFunction[1, 9][lna] * imPsi[lna]== 0,
		
		bFunction[2, 1][lna] * imDeltaPhi''[lna] + bFunction[2, 2][lna] * imDeltaPhi'[lna] + bFunction[2, 3][lna] * imDeltaPhi[lna] + bFunction[2, 4][lna] * imDeltaPsi''[lna] + bFunction[2, 5][lna] * imDeltaPsi'[lna] + bFunction[2, 6][lna] * imDeltaPsi[lna](* + bFunction[2, 7][lna] * imPsi''[lna]*) + bFunction[2, 8][lna] * imPsi'[lna] + bFunction[2, 9][lna] * imPsi[lna]== 0,
		
		bFunction[3, 1][lna] * imDeltaPhi''[lna] + bFunction[3, 2][lna] * imDeltaPhi'[lna] + bFunction[3, 3][lna] * imDeltaPhi[lna] + bFunction[3, 4][lna] * imDeltaPsi''[lna] + bFunction[3, 5][lna] * imDeltaPsi'[lna] + bFunction[3, 6][lna] * imDeltaPsi[lna](* + bFunction[3, 7][lna] * imPsi''[lna]*) + bFunction[3, 8][lna] * imPsi'[lna] + bFunction[3, 9][lna] * imPsi[lna]== 0
	};
	
	
	
	(* -------------------------------------------- *)
	
	(* Write down the initial conditions of deltaPhi(N) and deltaPsi(N) at the early times *)
	
	(* Here we omit the initial phase spectra. See Cosmology 28, p.49 and p.51. *)
	
	deltaPhiInitial = initialDeltaPhi;
	deltaPsiInitial = initialDeltaPsi;
	PsiInitial = initialPsi;
	
	
	
	(* -------------------------------------------- *)
	
	(* Calculate initial real and imaginary part *)
	
	(* Trivial when there's no initial phase *)
	
	reDeltaPhiInitial = Re[deltaPhiInitial];
	imDeltaPhiInitial = Im[deltaPhiInitial];
	
	reDeltaPsiInitial = Re[deltaPsiInitial];
	imDeltaPsiInitial = Im[deltaPsiInitial];
	
	rePsiInitial = Re[PsiInitial];
	imPsiInitial = Im[PsiInitial];
	
	
	
	(* -------------------------------------------- *)
	
	(* Write down the initial derivatives of deltaPhi and deltaPsi with respect to ln(a/ai) *)
	
	deltaPhiPrimeInitial = initialDeltaPhiPrime;
	
	deltaPsiPrimeInitial = initialDeltaPsiPrime;
	
	
	
	(* -------------------------------------------- *)
	
	(* Calculate initial derivatives of reDeltaPhi, imDeltaPhi, reDeltaPsi, and imDeltaPsi with respect to ln(a/ai) *)
	
	reDeltaPhiPrimeInitial = Re[deltaPhiPrimeInitial];
	imDeltaPhiPrimeInitial = Im[deltaPhiPrimeInitial];
	
	reDeltaPsiPrimeInitial = Re[deltaPsiPrimeInitial];
	imDeltaPsiPrimeInitial = Im[deltaPsiPrimeInitial];
	
	
	
	(* -------------------------------------------- *)
	
	(* Write down the initial conditions *)
	
	initialConditions =
	{
		reDeltaPhi[lnaBegin] == reDeltaPhiInitial,
		imDeltaPhi[lnaBegin] == imDeltaPhiInitial,
		reDeltaPhi'[lnaBegin] == reDeltaPhiPrimeInitial,
		imDeltaPhi'[lnaBegin] == imDeltaPhiPrimeInitial,
		
		reDeltaPsi[lnaBegin] == reDeltaPsiInitial,
		imDeltaPsi[lnaBegin] == imDeltaPsiInitial,
		reDeltaPsi'[lnaBegin] == reDeltaPsiPrimeInitial,
		imDeltaPsi'[lnaBegin] == imDeltaPsiPrimeInitial,
		
		rePsi[lnaBegin] == rePsiInitial,
		imPsi[lnaBegin] == imPsiInitial
	};
	
	
	
	(* -------------------------------------------- *)
	
	(* Record time for preparing *)
	
	timeEnd = DateList[];
	
	timePreparing = DateDifference[ timeBegin, timeEnd, "Second" ][[1]];
	
	(* Print["\nTime spent preparing (k = ", k * Mpc / PlanckLength, " \!\(\*SuperscriptBox[\(Mpc\), \(-1\)]\)) = ", timePreparing, " seconds."]; *)
	
	
	
	(* -------------------------------------------- *)
	
	(* Record time for solving *)
	
	timeBegin = DateList[];
	
	
	
	(* -------------------------------------------- *)
	
	(* NDSolve *)
	
	(* Mathematica default MaxSteps -> 10^4 *)
	
	sol =
		NDSolve[
			Join[differentialEquations, initialConditions],

			{reDeltaPhi, imDeltaPhi, reDeltaPsi, imDeltaPsi, rePsi, imPsi},
			{lna, lnaBegin, lnaEnd},

			MaxSteps -> 10^7
		];
	
	
	
	(* -------------------------------------------- *)
	
	(* Record time for solving *)
	
	timeEnd = DateList[];
	
	timeSolving = DateDifference[ timeBegin, timeEnd, "Second" ][[1]];
	
	(* Print["\nTime spent solving (k = ", k * Mpc / PlanckLength, " \!\(\*SuperscriptBox[\(Mpc\), \(-1\)]\)) = ", timeSolving, " seconds."]; *)
	
	
	
	(* -------------------------------------------- *)
	
	(* Retrieve the solutions *)
	
	reDeltaPhiFound = reDeltaPhi /. sol[[1]];
	imDeltaPhiFound = imDeltaPhi /. sol[[1]];
	reDeltaPsiFound = reDeltaPsi /. sol[[1]];
	imDeltaPsiFound = imDeltaPsi /. sol[[1]];
	rePsiFound = rePsi /. sol[[1]];
	imPsiFound = imPsi /. sol[[1]];
	
	lnaEndGot = reDeltaPhiFound[[ 1, 1, 2 ]];
	
	
	(* ///////////////////////////////////////////////
	(* -------------------------------------------- *)
	
	(* Plot the amplitude and the phase *)
	
	figDeltaPhiAmplitude =
		Plot[
			Abs[reDeltaPhiFound[lna] + I * imDeltaPhiFound[lna]],
			{lna, lnaBegin, lnaEndGot},
			
			PlotRange -> All,
			Axes -> False, Frame -> True, FrameTicks -> All,
			FrameLabel -> { "ln (a/\!\(\*SubscriptBox[\(a\), \(i\)]\))", "|\!\(\*SubscriptBox[\(\[Delta]\[Phi]\), \(k\)]\)\!\(\*SuperscriptBox[\(|\), \(2\)]\)/\!\(\*SubsuperscriptBox[\(m\), \(p\), \(-4\)]\)" },
			ImageSize -> 350
		];
	
	Print[
		"\n",
		"figDeltaPhiAmplitude (k = ", k * Mpc / PlanckLength, " \!\(\*SuperscriptBox[\(Mpc\), \(-1\)]\)) =\n",
		figDeltaPhiAmplitude
	];
	
	figDeltaPhiPhase =
		Plot[
			Arg[reDeltaPhiFound[lna] + I * imDeltaPhiFound[lna]],
			{lna, lnaBegin, lnaEndGot},
			
			PlotRange -> All,
			Axes -> False, Frame -> True, FrameTicks -> All,
			FrameLabel -> { "ln (a/\!\(\*SubscriptBox[\(a\), \(i\)]\))", "|\!\(\*SubscriptBox[\(\[Delta]\[Phi]\), \(k\)]\)\!\(\*SuperscriptBox[\(|\), \(2\)]\)/\!\(\*SubsuperscriptBox[\(m\), \(p\), \(-4\)]\)" },
			ImageSize -> 350
		];
	
	Print[
		"\n",
		"figDeltaPhiPhase (k = ", k * Mpc / PlanckLength, " \!\(\*SuperscriptBox[\(Mpc\), \(-1\)]\)) =\n",
		figDeltaPhiPhase
	];
	
	
	
	figDeltaPsiAmplitude =
		Plot[
			Abs[reDeltaPsiFound[lna] + I * imDeltaPsiFound[lna]],
			{lna, lnaBegin, lnaEndGot},
			
			PlotRange -> All,
			Axes -> False, Frame -> True, FrameTicks -> All,
			FrameLabel -> { "ln (a/\!\(\*SubscriptBox[\(a\), \(i\)]\))", "|\!\(\*SubscriptBox[\(\[Delta]\[Psi]\), \(k\)]\)\!\(\*SuperscriptBox[\(|\), \(2\)]\)/\!\(\*SubsuperscriptBox[\(m\), \(p\), \(-4\)]\)" },
			ImageSize -> 350
		];
	
	Print[
		"\n",
		"figDeltaPsiAmplitude (k = ", k * Mpc / PlanckLength, " \!\(\*SuperscriptBox[\(Mpc\), \(-1\)]\)) =\n",
		figDeltaPsiAmplitude
	];
	
	figDeltaPsiPhase =
		Plot[
			Arg[reDeltaPsiFound[lna] + I * imDeltaPsiFound[lna]],
			{lna, lnaBegin, lnaEndGot},
			
			PlotRange -> All,
			Axes -> False, Frame -> True, FrameTicks -> All,
			FrameLabel -> { "ln (a/\!\(\*SubscriptBox[\(a\), \(i\)]\))", "|\!\(\*SubscriptBox[\(\[Delta]\[Phi]\), \(k\)]\)\!\(\*SuperscriptBox[\(|\), \(2\)]\)/\!\(\*SubsuperscriptBox[\(m\), \(p\), \(-4\)]\)" },
			ImageSize -> 350
		];
	
	Print[
		"\n",
		"figDeltaPsiPhase (k = ", k * Mpc / PlanckLength, " \!\(\*SuperscriptBox[\(Mpc\), \(-1\)]\)) =\n",
		figDeltaPsiPhase
	];
	
	
	
	figPsiAmplitude =
		Plot[
			Abs[rePsiFound[lna] + I * imPsiFound[lna]],
			{lna, lnaBegin, lnaEndGot},
			
			PlotRange -> All,
			Axes -> False, Frame -> True, FrameTicks -> All,
			FrameLabel -> { "ln (a/\!\(\*SubscriptBox[\(a\), \(i\)]\))", "|\!\(\*SubscriptBox[\(\[Psi]\), \(k\)]\)\!\(\*SuperscriptBox[\(|\), \(2\)]\)/\!\(\*SubsuperscriptBox[\(m\), \(p\), \(-4\)]\)" },
			ImageSize -> 350
		];
	
	Print[
		"\n",
		"figPsiAmplitude (k = ", k * Mpc / PlanckLength, " \!\(\*SuperscriptBox[\(Mpc\), \(-1\)]\)) =\n",
		figPsiAmplitude
	];
	
	figPsiPhase =
		Plot[
			Arg[rePsiFound[lna] + I * imPsiFound[lna]],
			{lna, lnaBegin, lnaEndGot},
			
			PlotRange -> All,
			Axes -> False, Frame -> True, FrameTicks -> All,
			FrameLabel -> { "ln (a/\!\(\*SubscriptBox[\(a\), \(i\)]\))", "|\!\(\*SubscriptBox[\(\[Delta]\[Phi]\), \(k\)]\)\!\(\*SuperscriptBox[\(|\), \(2\)]\)/\!\(\*SubsuperscriptBox[\(m\), \(p\), \(-4\)]\)" },
			ImageSize -> 350
		];
	
	Print[
		"\n",
		"figPsiPhase (k = ", k * Mpc / PlanckLength, " \!\(\*SuperscriptBox[\(Mpc\), \(-1\)]\)) =\n",
		figPsiPhase
	];
	/////////////////////////////////////////////// *)
	
	
	(* -------------------------------------------- *)
	
	{reDeltaPhiFound, imDeltaPhiFound, reDeltaPsiFound, imDeltaPsiFound, rePsiFound, imPsiFound}
];



(* -------------------------------------------- *)
(* -------------------------------------------- *)

End[];



Protect[
	solvePerturbation
];

EndPackage[];
