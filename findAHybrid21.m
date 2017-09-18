(* :Title: Cosmology28/Numerical perturbation/findAHybrid21.m *)
(* :Notebook: Cosmology 23, p.109 *)

(* :Author: Yu-Hsiang Lin *)
(* :Date: Sun, Apr 6, 2014 *)

(* :Package Version: 2.1 *)
(* :Mathematica Version: 9.0 *)

(* :Context: findAHybrid21` *)

(*
:Requirements:


*)

(*
:Summary:

Find the dependences of the background fields on the scale factor in the inflation era.

Use Hybrid inflation model (Copeland notation).

Note that all quantities in this program are in
	Planck units.
*)

(*
:Keywords:


*)

(*
:Discussion:

Version 2.0:
	Change the formulation to use ln(a/ai) as the variable.

Version 2.1:
	Change MaxSteps in NDSolve[].
*)



(* -------------------------------------------- *)
(* -------------------------------------------- *)

BeginPackage["findAHybrid21`"];



(* -------------------------------------------- *)

(* Functions to be exported *)

{
	findAHybrid
};



(* -------------------------------------------- *)
(* -------------------------------------------- *)

Begin["`Private`"];



(* -------------------------------------------- *)
(* Solve inflation era. *)
(* -------------------------------------------- *)

(* NDSolve *)

(* We don't need to specify a and t anymore. Only the final ln(a) is needed. *)

findAHybrid[potential_, initialPhiBar_, initialDPhiBarDt_, initialPsiBar_, initialDPsiBarDt_, finallna_]:=
Module[
	{
		lna,
		X,
		differentialEquations, initialConditions, sol,
		phiBar, psiBar,
		phiBarFound, psiBarFound
	},
	
	Off[NDSolve::mxst];
	
	(* The kinetic energy term expressed in phiBar'[lna] and psiBar'[lna]. *)
	
	X[lnaVar_] =
		1/2 * (
			8 * Pi * potential[phiBar[lnaVar], psiBar[lnaVar]] * ( phiBar'[lnaVar]^2 + psiBar'[lnaVar]^2 ) /
			(  3 - 4 * Pi * ( phiBar'[lnaVar]^2 + psiBar'[lnaVar]^2 )  )
		);
	
	differentialEquations =
	{
		phiBar''[lna] +
		3 * potential[phiBar[lna], psiBar[lna]] / (  X[lna] + potential[phiBar[lna], psiBar[lna]]  ) * phiBar'[lna] +
		3 / (  8 * Pi * ( X[lna] + potential[phiBar[lna], psiBar[lna]] )  ) * D[ potential[phiBar[lna], psiBar[lna]], phiBar[lna] ] == 0,
		
		psiBar''[lna] +
		3 * potential[phiBar[lna], psiBar[lna]] / (  X[lna] + potential[phiBar[lna], psiBar[lna]]  ) * psiBar'[lna] +
		3 / (  8 * Pi * ( X[lna] + potential[phiBar[lna], psiBar[lna]] )  ) * D[ potential[phiBar[lna], psiBar[lna]], psiBar[lna] ] == 0
	};



	initialConditions =
	{
		phiBar[0] == initialPhiBar,
		phiBar'[0] == ( 3 / ( 8 * Pi ) )^(1/2) * initialDPhiBarDt / ( initialDPhiBarDt^2 / 2 + initialDPsiBarDt^2 / 2 + potential[initialPhiBar, initialPsiBar] )^(1/2),
		
		psiBar[0] == initialPsiBar,
		psiBar'[0] == ( 3 / ( 8 * Pi ) )^(1/2) * initialDPsiBarDt / ( initialDPhiBarDt^2 / 2 + initialDPsiBarDt^2 / 2 + potential[initialPhiBar, initialPsiBar] )^(1/2)
	};

	sol =
		NDSolve[
			Join[differentialEquations, initialConditions],

			{phiBar, psiBar},
			{lna, 0, finallna},

			MaxSteps -> 10^6   (* 10^7 *)
		];

	On[NDSolve::mxst];

	(* -------------------------------------------- *)

	phiBarFound = phiBar /. sol[[1]];
	psiBarFound = psiBar /. sol[[1]];



	(* -------------------------------------------- *)

	{phiBarFound, psiBarFound}
];



(* -------------------------------------------- *)
(* -------------------------------------------- *)

End[];



Protect[
	findAInflation
];

EndPackage[];
