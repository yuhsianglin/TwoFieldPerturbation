(* :Title: Cosmology36/NumericalPerturbation/Main560.m *)
(* :Notebook: Cosmology 36, p.31 *)

(* :Author: Yu-Hsiang Lin *)
(* :Date: Thu, Apr 16, 2015 *)

(* :Script Version: 55.0 *)
(* :Mathematica Version: 9.0 *)

(* :Context: Global` *)

(*
:Requirements:

Cosmology30/NumericalPerturbation/ZeroOrderPartXXX.m
Cosmology30/NumericalPerturbation/FirstOrderPartXXX.m

Cosmology29/NumericalPerturbation/findAHybrid21.m
Cosmology27/CambrianExplosion/FUN30.m

Cosmology30/NumericalPerturbation/solveHybridPerturbation203.m
*)

(*
:Summary:

The main program of the parameter scanning.
*)

(*
:Keywords:


*)

(*
:Discussion:

Version 32.0:
	Automate.
	This is the first version of the Main-series.

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



(* ------------------------------------------- *)
(* Load packages *)
(* ------------------------------------------- *)

Needs["findAHybrid21`"];
Needs["FUN30`"];

Needs["solveHybridPerturbation203`"];



(* ------------------------------------------- *)
(* Scan lambdaPrime and then initialPhiBar, with initialPsiBar fixed *)
(* ------------------------------------------- *)

(* The folder where our source codes locate *)

calculationFolder =
	"/Users/zachlin/Dropbox/Mathematica201501/Cosmology36/NumericalPerturbation";



(* The folder where we store the results *)

outputFolder =
	CreateDirectory["/Users/zachlin/Dropbox/Mathematica201501/Cosmology36/NumericalPerturbation/Output560"];



(* ------------------------------------------- *)

(* Specify the range of parameter scan *)

lambdaPrimeList =
	{1.*10^-9};

(*	{1.*10^-10, 1.*10^-9, 3.*10^-9};*)

mList =
	Table[
		mValue,
		
		{mValue, 1.22*10^-6, 1.22*10^-6, 0.002*10^-6}
	];

initialPsiBarList =
	Table[
		initialPsiBarValue,
		
		{initialPsiBarValue, 1.64, 1.64, 0.1}
	];



Print[
	"\n",
	"lambdaPrimeList = ", lambdaPrimeList
];

Print[
	"\n",
	"mList = ", mList
];

Print[
	"\n",
	"initialPsiBarList = ", initialPsiBarList
];



(* ------------------------------------------- *)

(* Calculate *)

SetDirectory[calculationFolder];

Do[
(
	lambdaPrime = lambdaPrimeList[[lambdaPrimeIndex]];
	m = mList[[mIndex]];
	initialPsiBar = initialPsiBarList[[initialPsiBarIndex]];
	
	<< ZeroOrderPart560.m;
	<< FirstOrderPart560.m;
),
	
	{lambdaPrimeIndex, 1, Length[lambdaPrimeList]},
	{mIndex, 1, Length[mList]},
	{initialPsiBarIndex, 1, Length[initialPsiBarList]}
];

ResetDirectory[];
