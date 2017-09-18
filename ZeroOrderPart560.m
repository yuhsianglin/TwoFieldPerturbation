(* :Title: Cosmology36/NumericalPerturbation/ZeroOrderPart560.m *)
(* :Notebook: Cosmology 36, p.31 *)

(* :Author: Yu-Hsiang Lin *)
(* :Date: Thu, Apr 16, 2015 *)

(* :Script Version: 55.0 *)
(* :Mathematica Version: 9.0 *)

(* :Context: Global` *)

(*
:Requirements:

Cosmology29/NumericalPerturbation/findAHybrid21.m
Cosmology27/CambrianExplosion/FUN30.m
*)

(*
:Summary:

Find phiBar(lna), psiBar(lna).

Note that all quantities in this program are in
	Planck units.

Here the scale factor a should be understood as
	atilde = a / ai.
*)

(*
:Keywords:


*)

(*
:Discussion:

Version 1.0:
	

Version 2.0:
	We change the formulation to use ln(a/ai) as our variable, and no longer directly solve for a(t) anymore.

Version 3.0:
	Update to solveHybridPerturbation20.m.

Version 4.0:
	Update to solveHybridPerturbation30.m.

Version 5.0:
	Update to solveHybridPerturbation40.m.

Version 5.1:
	Update to solveHybridPerturbation41.m.

Version 6.0:
	Minor changes:
		Add phiBarFound and psiBarFound.
		Make plots of phiBarFound and psiBarFound.
		Print e-Fold.

Version 7.0:
	Minor changes:
		In phase diagram, \psi is in x axis and \phi is in y axis.
		Migrate to Mathematica 9.

Version 8.0:
	Same as 7.0.
	Delta-N formalism.

Version 8.1:
	Alter lambda, lambdaPrime, phi_i, psi_i.

Version 9.0:
	Using "dynamical H".
	Trying out different parameters.

Version 9.2:
	Checking the horizon exit time.
	Use FUN30.m.

Version 10.0:
	Change initial conditions in perturbation.

Version 11.0:
	Check the process for horizon exit.
	In the .nb file, I fix the connection between transient and slow-roll periods.

Version 12.0:
	Check early time approximation.

Version 13.0:
	Solve again.

Version 14.0:
	Improve the solution.
	Use "findAHybrid21`".

Version 14.1:
	Use "findAHybrid20`".

Version 15.0:
	Bridge psi-dominated and phi-dominated eras.

Version 15.1:
	Bridge psi-dominated and phi-dominated eras.
	Use "findAHybrid21`".

Version 15.2:
	Bridge psi-dominated and phi-dominated eras.

Version 16.0:
	Calculate the perturbation.
	Now combine the single-field inflation regime (stage 3) into the script.
	Verify that the well-within-horizon approximation works well.

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
	Change psiBar3[lna] into psiBarFound3[lna] so that the naming is more consistent.

Version 21.0:
	Solve for the spectrum.
	Use Piecewise[].

Version 22.0:
	Small interaction strength.
	Increase the field values.

Version 22.1:
	Tune.

Version 25.0:
	Try out a specific parameter set.

Version 26.0:
	Try to move the effect to larger scales.

Version 26.1:
	Try larger inflaton mass.

Version 27.X:
	Scan initialPsiBar.

Version 28.X:
	Scan initialPhiBar.

Version 29.X:
	Scan lambdaPrime.

Version 30.0:
	Find how far it can go for super-horizon modes.

Version 31.0:
	Check whether the choice of an earlier N to begin solving the perturbation changes the outcome.

Version 32.0:
	Automate.

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

Needs["findAHybrid21`"];

Needs["FUN30`"];
////////////////////////////////////////////// *)



(* ------------------------------------------- *)
(* Solve phiBar, psiBar: Stage 1, 2, and 3 *)
(* ------------------------------------------- *)

(* Potential: Copeland notation *)

potential[phi_, psi_] = 
	lambdaPrime/2 * phi^2 * psi^2 + 1/2 * m^2 * phi^2;



(* Model parameters *)

(* lambda = 1.*10^-9; *)
(* lambdaPrime = 1.*10^-9; *)   (* Now it's the parameter to scan *)
(* m = 1.3*10^-6; *)   (* Now it's the parameter to scan *)
(* M = 1.*10^-3; *)



(* ------------------------------------------- *)

(* Initial conditions *)

initialPhiBar = 3.65;
initialDPhiBarDt = -1.*10^-8;

(*initialPsiBar = 1.6;*)
initialDPsiBarDt = -1.*10^-8;



(* Final ln(a) *)

finallna = 200.;



(* ------------------------------------------- *)

(* Solve *)

(* Note that the output of findAHybrid[] is {phiBar(ln a), psiBar(ln a)} *)

solution =
	findAHybrid[potential, initialPhiBar, initialDPhiBarDt, initialPsiBar, initialDPsiBarDt, finallna];



(* ------------------------------------------- *)

(* Retrieve solutions *)

phiBarFound = solution[[1]];
psiBarFound = solution[[2]];



(* ------------------------------------------- *)
(* Find the e-fold number when the energy density drops to the reheating energy density *)
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



n34Trial = 70;

n34 = n /. FindRoot[3 * H[n]^2 / (8*Pi) == 2.67523 * 10^(-13), {n, n34Trial}];
(* 2.67523 * 10^(-13) *)


Print[
	"\n",
	"--------------------------------------------------------", "\n",
	"\n",
	"m = ", m, "\n",
	"lambdaPrime = ", lambdaPrime, "\n",
	"initialPhiBar = ", initialPhiBar, "\n",
	"initialPsiBar = ", initialPsiBar, "\n"
];



Print[
	"\n",
	"The total number of e-folds is given by\n",
	"n34 = ", n34
];



(* ------------------------------------------- *)
(* Plot phiBarFound and psiBarFound *)
(* ------------------------------------------- *)

figPhiBar = 
	Plot[phiBarFound[lna], {lna, 0, n34}, PlotRange -> All, Axes -> False, Frame -> True, FrameTicks -> All, FrameLabel -> {"ln(a/ai)", "\[Phi]/\!\(\*SubscriptBox[\"m\", \"p\"]\)"}, ImageSize -> 350];

Print[
	"\n",
	"figPhiBar = \n",
	figPhiBar
];



figPsiBar = 
	Plot[psiBarFound[lna], {lna, 0, n34}, PlotRange -> All, Axes -> False, Frame -> True, FrameTicks -> All, FrameLabel -> {"ln(a/ai)", "\[Psi]/\!\(\*SubscriptBox[\"m\", \"p\"]\)"}, ImageSize -> 350];

Print[
	"\n",
	"figPsiBar = \n",
	figPsiBar
];
