(* :Title: Cosmology27/CambrianExplosion/FUN30.m *)
(* :Notebook: None *)

(* :Author: Yu-Hsiang Lin *)
(* :Date: Sun, Feb 16, 2014 *)

(* :Package Version: 3.0 *)
(* :Mathematica Version: 9.0 *)

(* :Context: FUN30` *)

(*
:Requirements:


*)

(*
:Summary:

	Frequently Used Numbers.
*)

(*
:Keywords:


*)

(*
:Discussion:

Version 2.0:
	Add Planck mass.

Version 3.0:
	Add vacuum permittivity, electron mass, electronRadius, barn, and erg.
*)



(* -------------------------------------------- *)
(* -------------------------------------------- *)

BeginPackage["FUN30`"];


(*
(* -------------------------------------------- *)

(* Functions to be exported *)

{
	
};
*)


(* -------------------------------------------- *)

(* Basic constants *)

c = 2.99792458*10^8;
hbar = 1.054572*10^-34;
kB = 1.38065*10^-23;
G = 6.67428*10^-11;
electronCharge = 1.60218*10^-19;

Mpc = 30.857*10^15*10^6;
year = 86400*365.25;

PlanckMass = ( hbar * c / G )^(1/2);
PlanckMassInGeV = PlanckMass * c^2 / ( 10^9 * electronCharge );

PlanckTime = ( hbar * G / c^5)^(1/2);

PlanckLength = ( hbar * G / c^3 )^(1/2);

epsilon0 = 8.8541878176 * 10^-12;
electronMass = 9.109382 * 10^-31;
electronRadius = electronCharge^2 / ( 4 * Pi * epsilon0 * electronMass * c^2 );

barn = 10^-28;
erg = 10^-7;


(*
(* -------------------------------------------- *)
(* -------------------------------------------- *)

Begin["`Private`"];



(* -------------------------------------------- *)
(* -------------------------------------------- *)

End[];
*)


Protect[
	c, hbar, kB, G, electronCharge,
	Mpc, year,
	PlanckMass, PlanckMassInGeV, PlanckTime, PlanckLength,
	epsilon0, electronMass, electronRadius,
	barn, erg
];

EndPackage[];
