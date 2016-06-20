(* ::Package:: *)

$Assumptions = {e >= m, m >= 0, 0 <= \[Theta] < \[Pi], 0 <= \[Phi] < 2\[Pi]};
MergeMatrix[m_] := Join @@ (MapThread[Join, #]& /@ m)
I2 = IdentityMatrix[2];
I4 = IdentityMatrix[4];
O2 = 0*I2;
O4 = 0*I4;
XtoPhi[exp_] := exp //. {x -> Sqrt[e-m] Sqrt[e+m]Sin[\[Theta]]Cos[\[Phi]], y -> Sqrt[e-m] Sqrt[e+m]Sin[\[Theta]]Sin[\[Phi]], z -> Sqrt[e-m] Sqrt[e+m]Cos[\[Theta]]} // Simplify


(* ::Subsubsection:: *)
(*Two component spinor*)


\[Sigma][0]   = I2;
\[Sigma][i_] := PauliMatrix[i]
J[i_] := \[Sigma][i]/2
\[Gamma][0]        = MergeMatrix[{{I2, O2  }, { O2  , -I2}}];
\[Gamma][i:1|2|3] := MergeMatrix[{{O2, \[Sigma][i]}, {-\[Sigma][i],  O2}}];
\[Gamma][5]       := I \[Gamma][0].\[Gamma][1].\[Gamma][2].\[Gamma][3]

S[i_] := (1/2) MergeMatrix[{{\[Sigma][i], O2},{O2, \[Sigma][i]}}]
\[CapitalLambda][J_, \[Theta]_, \[Phi]_]  := {Sin[\[Theta]]Cos[\[Phi]], Sin[\[Theta]]Sin[\[Phi]], Cos[\[Theta]]}.{J[1], J[2], J[3]};
U2[\[Theta]_, \[Phi]_] := MatrixExp[-I \[Phi] J[3]].MatrixExp[-I \[Theta] J[2]].MatrixExp[I \[Phi] J[3]]
U4[\[Theta]_, \[Phi]_] := MatrixExp[-I \[Phi] S[3]].MatrixExp[-I \[Theta] S[2]].MatrixExp[I \[Phi] S[3]]
\[Chi][z,  1/2] = {{1}, {0}};
\[Chi][z, -1/2] = {{0}, {1}};
(*1.9*)  \[Chi][ p, \[Lambda]_] := U2[\[Theta], \[Phi]].\[Chi][z, \[Lambda]]
(*1.19*) \[Chi][-z, \[Lambda]_] := (-1)^(1/2-\[Lambda])MatrixExp[-I \[Pi] J[2]].\[Chi][z, \[Lambda]]
(*1.21*) \[Chi][-p, \[Lambda]_] := U2[\[Theta], \[Phi]].\[Chi][-z, \[Lambda]]


(* ::Text:: *)
(*Note that (1.19) and (1.21) is based on Jacob-Wick convention, i.e., so that  \[Chi](-z, \[Lambda]) = \[Chi](z, -\[Lambda]).*)
(*(Cf. Eq. 13 of Jacob and Wick, Annal. Phys. 7 404-428)*)
(*Then, obviously, \[Chi](-p, \[Lambda]) defined by rotation from \[Chi](-z, \[Lambda]) satisfies \[Chi](-p, \[Lambda]) = \[Chi](p, -\[Lambda]).*)


(*1.20*)
\[Chi][-z, #] == \[Chi][z, -#] &/@ {1/2, -1/2}
\[Chi][-p, #] == \[Chi][p, -#] &/@ {1/2, -1/2}


(* ::Text:: *)
(*The phase difference between Jacob-Wick convention and a naive replacement of (\[Theta], \[Phi])->(\[Pi]-\[Theta], \[Phi]+\[Pi]):*)


(((2# Exp[-2 I # \[Phi]]\[Chi][p, #]) /. {\[Theta] -> \[Pi]-\[Theta], \[Phi] -> \[Pi]+\[Phi]}) == \[Chi][-p, #]) &/@ {1/2, -1/2} // FullSimplify


(* ::Text:: *)
(*Note that this definition is consistent with helicity:*)


(\[CapitalLambda][J, 0,   \[Phi]  ].\[Chi][ z, #] == # \[Chi][ z, #]) &/@ {1/2, -1/2} // FullSimplify
(\[CapitalLambda][J, \[Pi],   \[Phi]  ].\[Chi][-z, #] == # \[Chi][-z, #]) &/@ {1/2, -1/2} // FullSimplify
(\[CapitalLambda][J, \[Theta],   \[Phi]  ].\[Chi][ p, #] == # \[Chi][ p, #]) &/@ {1/2, -1/2} // FullSimplify
(\[CapitalLambda][J, \[Pi]-\[Theta], \[Phi]+\[Pi]].\[Chi][-p, #] == # \[Chi][-p, #]) &/@ {1/2, -1/2} // FullSimplify
-\[CapitalLambda][J, \[Theta], \[Phi]] == \[CapitalLambda][J, \[Pi]-\[Theta], \[Phi]+\[Pi]]


(* ::Text:: *)
(*Several validations:*)


(*1.38*)
\[Chi][p, 1/2]
\[Chi][p, -1/2]
(*1.10*)
(J[3].\[Chi][z, #] == # \[Chi][z, #]) &/@ {1/2, -1/2} // FullSimplify
(*1.39*)
\[Chi][p, #] == MatrixExp[-I \[Theta] (Cos[\[Phi]]J[2]-Sin[\[Phi]]J[1])].\[Chi][z, #] &/@ {1/2, -1/2} // FullSimplify
\[Chi][-p, #] == MatrixExp[-I \[Theta] (Cos[\[Phi]]J[2]-Sin[\[Phi]]J[1])].\[Chi][-z, #] &/@ {1/2, -1/2} // FullSimplify


(* ::Text:: *)
(*Eq.(1.40) holds only for the first particle. The second particle in Jacob-Wick convention requires extra factor (-1).*)


(*1.40*)
Correction140[P:z|p] := 1
Correction140[P:-z|-p] := -1
Do[Correction140[P] \[Chi][P, -\[Lambda]] == -2 \[Lambda] I \[Sigma][2].Conjugate[\[Chi][P, \[Lambda]]] // FullSimplify // Print, {P, {z, p, -z, -p}}, {\[Lambda], {1/2, -1/2}}]


(* ::Subsubsection:: *)
(*Four component spinor*)


(* ::Text:: *)
(*The u-spinor is defined by (1.34), which is equivalent to (1.41a) as we see below.*)
(*The v-spinor is by (1.36), but (1.41b) is not valid for the second particle because of error in (1.40).*)


bar[spinor_] := Conjugate[Transpose[spinor]].\[Gamma][0] // FullSimplify

(*1.41a*) u[p_, \[Lambda]_] := MergeMatrix[{{Sqrt[e+m] \[Chi][p, \[Lambda]]}, { 2\[Lambda] Sqrt[e-m] \[Chi][p, \[Lambda]]}}]
(*1.36*)  v[p_, \[Lambda]_] := I \[Gamma][0].\[Gamma][2].Transpose[bar[u[p, \[Lambda]]]]


Module[{\[Eta], p\[Sigma], mom = Sqrt[e-m] Sqrt[e+m]},
  \[Eta][p_, \[Lambda]_] := -I \[Sigma][2].Conjugate[\[Chi][p, \[Lambda]]] // Simplify;
  p\[Sigma][z] := \[Sigma][3]mom;
  p\[Sigma][-z] := -\[Sigma][3]mom;
  p\[Sigma][p]  :=  {Sin[\[Theta]]Cos[\[Phi]], Sin[\[Theta]]Sin[\[Phi]], Cos[\[Theta]]}.{\[Sigma][1], \[Sigma][2], \[Sigma][3]}mom;
  p\[Sigma][-p] := -{Sin[\[Theta]]Cos[\[Phi]], Sin[\[Theta]]Sin[\[Phi]], Cos[\[Theta]]}.{\[Sigma][1], \[Sigma][2], \[Sigma][3]}mom;
  (*1.34*)
  Print["* Peskin u(p,\[Lambda])"];
  Do[u[P, \[Lambda]] == MergeMatrix[{{Sqrt[e+m] \[Chi][P, \[Lambda]]}, {1/Sqrt[e+m] p\[Sigma][P].\[Chi][P, \[Lambda]]}}] // Simplify // Print, {P, {z, p, -z, -p}}, {\[Lambda], {1/2, -1/2}}];
  Print["* Peskin v(p,\[Lambda])"];
  Do[v[P, \[Lambda]] == MergeMatrix[{{-(1/Sqrt[e+m])p\[Sigma][P].\[Eta][P, \[Lambda]]}, {-Sqrt[e+m] \[Eta][P, \[Lambda]]}}] // Simplify // Print, {P, {z, p, -z, -p}}, {\[Lambda], {1/2, -1/2}}]
  Print["* 9405376 (1.41)"];
  Do[Correction140[P]v[P, \[Lambda]] == MergeMatrix[{{Sqrt[e-m] \[Chi][P, -\[Lambda]]}, {-2\[Lambda] Sqrt[e+m] \[Chi][P, -\[Lambda]]}}] // Print, {P, {z, p, -z, -p}}, {\[Lambda], {1/2, -1/2}}]
]


(* ::Text:: *)
(*This is consistent with helicity, and satisfies Jacob-Wick condition.*)
(*(Note that v(\[Lambda]) has a helicity -\[Lambda].)*)


-\[CapitalLambda][S, \[Theta], \[Phi]] == \[CapitalLambda][S, \[Pi]-\[Theta], \[Phi]+\[Pi]]
(\[CapitalLambda][S, 0,   \[Phi]  ].u[ z, #] == # u[ z, #]) &/@ {1/2, -1/2} // FullSimplify
(\[CapitalLambda][S, \[Pi],   \[Phi]  ].u[-z, #] == # u[-z, #]) &/@ {1/2, -1/2} // FullSimplify
(\[CapitalLambda][S, \[Theta],   \[Phi]  ].u[ p, #] == # u[ p, #]) &/@ {1/2, -1/2} // FullSimplify
(\[CapitalLambda][S, \[Pi]-\[Theta], \[Phi]+\[Pi]].u[-p, #] == # u[-p, #]) &/@ {1/2, -1/2} // FullSimplify
(\[CapitalLambda][S, 0,   \[Phi]  ].v[ z, #] == -# v[ z, #]) &/@ {1/2, -1/2} // FullSimplify
(\[CapitalLambda][S, \[Pi],   \[Phi]  ].v[-z, #] == -# v[-z, #]) &/@ {1/2, -1/2} // FullSimplify
(\[CapitalLambda][S, \[Theta],   \[Phi]  ].v[ p, #] == -# v[ p, #]) &/@ {1/2, -1/2} // FullSimplify
(\[CapitalLambda][S, \[Pi]-\[Theta], \[Phi]+\[Pi]].v[-p, #] == -# v[-p, #]) &/@ {1/2, -1/2} // FullSimplify

(u[z, #] == u[p, #] /. \[Theta] -> 0) &/@ {1/2, -1/2}
(v[z, #] == v[p, #] /. \[Theta] -> 0) &/@ {1/2, -1/2}
(u[-z, #] == u[-p, #] /. \[Theta] -> 0) &/@ {1/2, -1/2}
(v[-z, #] == v[-p, #] /. \[Theta] -> 0) &/@ {1/2, -1/2}
u[p, -1/2] == u[-p, 1/2] /. {e->m}
v[z, -1/2] == v[-z, 1/2] /. {e->m}
v[p, -1/2] == v[-p, 1/2] /. {e->m}


(* ::Text:: *)
(*Formulae (1.42) is correct, while (1.43) and (1.44b) need modification for the second particle.*)


Correction143[P: z | p] := 1
Correction143[P: -z|-p] := -1
Correction144b[_] := -1
(*1.42*)
\[Chi][-z, -#] == \[Chi][z, #] &/@ {1/2, -1/2}
\[Chi][-p, -#] == \[Chi][p, #] &/@ {1/2, -1/2}
(*1.43*)
Do[Correction143[P] v[P, \[Lambda]] == -2 \[Lambda] \[Gamma][5].u[P, -\[Lambda]] // Simplify // Print, {P, {z, p, -z, -p}}, {\[Lambda], {1/2, -1/2}}]
(*1.44*)
Do[u[-P, -\[Lambda]] == \[Gamma][0].u[P, \[Lambda]] // Simplify // Print, {P, {z, p, -z, -p}}, {\[Lambda], {1/2, -1/2}}]
Do[Correction144b[P] v[-P, -\[Lambda]] == \[Gamma][0].v[P, \[Lambda]] // Simplify // Print, {P, {z, p, -z, -p}}, {\[Lambda], {1/2, -1/2}}]


(* ::Subsubsection:: *)
(*Helicity sum *)


s[ z, \[Lambda]_] := (2\[Lambda]/m){{Sqrt[e^2-m^2]}, {0}, {0},  {e}}
s[ p, \[Lambda]_] := (2\[Lambda]/m){{Sqrt[e^2-m^2]}, {e Sin[\[Theta]]Cos[\[Phi]]}, {e Sin[\[Theta]]Sin[\[Phi]]}, {e Cos[\[Theta]]}}
s[-z, \[Lambda]_] := (2\[Lambda]/m){{Sqrt[e^2-m^2]}, {0}, {0}, {-e}}
s[-p, \[Lambda]_] := (2\[Lambda]/m){{Sqrt[e^2-m^2]}, {-e Sin[\[Theta]]Cos[\[Phi]]}, {-e Sin[\[Theta]]Sin[\[Phi]]}, {-e Cos[\[Theta]]}}
FourDot[a_List, b_List] := Transpose[a].DiagonalMatrix[{1, -1, -1, -1}].b;


(*1.49*)
FourDot[s[z, \[Lambda]], {{e}, {0}, {0}, {Sqrt[e^2-m^2]}}] == 0 // Simplify
FourDot[s[p, \[Lambda]], {{e}, {Sqrt[e^2-m^2]Sin[\[Theta]]Cos[\[Phi]]}, {Sqrt[e^2-m^2]Sin[\[Theta]]Sin[\[Phi]]}, {Sqrt[e^2-m^2]Cos[\[Theta]]}}] == 0 // Simplify
FourDot[s[z, #], s[z, #]] == -1 &/@ {1/2, -1/2} // Simplify
FourDot[s[p, #], s[p, #]] == -1 &/@ {1/2, -1/2} // Simplify


SlashP[z]  = FourDot[{{e}, {0}, {0},  { Sqrt[e^2-m^2]}}, {{\[Gamma][0]}, {\[Gamma][1]}, {\[Gamma][2]}, {\[Gamma][3]}}][[1,1]];
SlashP[-z] = FourDot[{{e}, {0}, {0},  {-Sqrt[e^2-m^2]}}, {{\[Gamma][0]}, {\[Gamma][1]}, {\[Gamma][2]}, {\[Gamma][3]}}][[1,1]];
SlashP[p]  = FourDot[{{e}, {Sqrt[e^2-m^2] Sin[\[Theta]]Cos[\[Phi]]}, {Sqrt[e^2-m^2] Sin[\[Theta]]Sin[\[Phi]]}, {Sqrt[e^2-m^2] Cos[\[Theta]]}}, {{\[Gamma][0]}, {\[Gamma][1]}, {\[Gamma][2]}, {\[Gamma][3]}}][[1,1]];
SlashP[-p] = FourDot[{{e}, {-Sqrt[e^2-m^2] Sin[\[Theta]]Cos[\[Phi]]}, {-Sqrt[e^2-m^2] Sin[\[Theta]]Sin[\[Phi]]}, {-Sqrt[e^2-m^2] Cos[\[Theta]]}}, {{\[Gamma][0]}, {\[Gamma][1]}, {\[Gamma][2]}, {\[Gamma][3]}}][[1,1]];
SlashS[P_, \[Lambda]_] := FourDot[s[P, \[Lambda]], {{\[Gamma][0]}, {\[Gamma][1]}, {\[Gamma][2]}, {\[Gamma][3]}}][[1,1]]


(*1.50*)
Do[(SlashP[P].u[P, \[Lambda]] ==  m u[P, \[Lambda]]) // Simplify // Print, {P, {z, p, -z, -p}}, {\[Lambda], {1/2, -1/2}}]
Do[(SlashP[P].v[P, \[Lambda]] == -m v[P, \[Lambda]]) // Simplify // Print, {P, {z, p, -z, -p}}, {\[Lambda], {1/2, -1/2}}]
Do[(\[Gamma][5].SlashS[P, \[Lambda]].u[P, \[Lambda]] == u[P, \[Lambda]]) // Simplify // Print, {P, {z, p, -z, -p}}, {\[Lambda], {1/2, -1/2}}]
Do[(\[Gamma][5].SlashS[P, \[Lambda]].v[P, \[Lambda]] == v[P, \[Lambda]]) // Simplify // Print, {P, {z, p, -z, -p}}, {\[Lambda], {1/2, -1/2}}]


(*1.51*)
Do[(\[Gamma][5].SlashS[P, \[Lambda]].u[P, \[Lambda]] == u[P, \[Lambda]]) // Simplify // Print, {P, {z, p, -z, -p}}, {\[Lambda], {1/2, -1/2}}]
Do[(\[Gamma][5].SlashS[P, \[Lambda]].v[P, \[Lambda]] == v[P, \[Lambda]]) // Simplify // Print, {P, {z, p, -z, -p}}, {\[Lambda], {1/2, -1/2}}]


Do[(u[P, \[Lambda]].bar[u[P, \[Lambda]]] == (1/2)(I4 + \[Gamma][5].SlashS[P, \[Lambda]]).(SlashP[P] + m I4)) // FullSimplify // Print, {P, {z, p, -z, -p}}, {\[Lambda], {1/2, -1/2}}]
Do[(v[P, \[Lambda]].bar[v[P, \[Lambda]]] == (1/2)(I4 + \[Gamma][5].SlashS[P, \[Lambda]]).(SlashP[P] - m I4)) // FullSimplify // Print, {P, {z, p, -z, -p}}, {\[Lambda], {1/2, -1/2}}]



