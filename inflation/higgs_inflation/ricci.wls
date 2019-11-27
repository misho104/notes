#!/usr/bin/env wolframscript
(* ::Package:: *)

SetDirectory[NotebookDirectory[]];
<<WeylAlgebra`


Delete[MakeBoxesNT$deriv];
MakeBoxesNT$deriv[f_, c_, Field[n_], any___] := MakeBoxesNT$deriv[f, Style[#, Blue]&, n, any]
MakeBoxesNT$deriv[f_, c_, Field[Dx[n_]], UI[\[Mu]_, "lorentz"], any___] := Row[{Superscript["\[PartialD]", \[Mu]], MakeBoxesNT$deriv[f, c, Field[n], any]}]
MakeBoxesNT$deriv[f_, c_, Field[Dx[n_]], LI[\[Mu]_, "lorentz"], any___] := Row[{Subscript["\[PartialD]", \[Mu]], MakeBoxesNT$deriv[f, c, Field[n], any]}]
MakeBoxesNT$deriv[f_, c_, OverBar[Field[Dx[n_]]], UI[\[Mu]_, "lorentz"], any___] := Row[{Superscript["\[PartialD]", \[Mu]], MakeBoxesNT$deriv[f, c, OverBar[Field[n]], any]}]
MakeBoxesNT$deriv[f_, c_, OverBar[Field[Dx[n_]]], LI[\[Mu]_, "lorentz"], any___] := Row[{Subscript["\[PartialD]", \[Mu]], MakeBoxesNT$deriv[f, c, OverBar[Field[n]], any]}]
MakeBoxesNT$deriv[f_, c_, Field[Conj[n_[i_]]], any___] := MakeBoxesNT$deriv[f, c, Field[Subsuperscript[n,i,"*"]], any]
MakeBoxesNT$deriv[f_, c_, Field[Conj[n_]], any___] := MakeBoxesNT$deriv[f, c, Field[Superscript[n,"*"]], any]
MakeBoxesNT$deriv[f_, c_, Field[n_[i_]], any___] := MakeBoxesNT$deriv[f, c, Subscript[n,i], any]
MakeBoxesNT$deriv[f_, c_, Field[n_], any___] := MakeBoxesNT$deriv[f, c, n, any]
MakeBoxesNT$deriv[f_, c_, OverBar[Field[n_[i_]]], any___] := MakeBoxesNT$deriv[f, c, OverBar[Subscript[n,i]], any]
MakeBoxesNT$deriv[f_, c_, OverBar[Field[n_]], any___] := MakeBoxesNT$deriv[f, c, OverBar[n], any]

MakeBoxesNT$deriv[f_, c_, Field[Dx2[n_]], any___] := Row[{Superscript["\[PartialD]", 2], MakeBoxesNT$deriv[f, c, Field[n], any]}]
MakeBoxesNT$deriv[exp___] := WeylAlgebra`Private`MakeBoxesNT$[exp]
WeylAlgebra`MakeBoxesNT = MakeBoxesNT$deriv;

Deriv /: MakeBoxes[obj: Deriv["x"[UI[i_,"lorentz"]], exp_], f:StandardForm|TraditionalForm] := Row[{Superscript["\[PartialD]", i], exp}] // ToBoxes // InterpretationBox[#,obj] &
Deriv /: MakeBoxes[obj: Deriv["x"[LI[i_,"lorentz"]], exp_], f:StandardForm|TraditionalForm] := Row[{Subscript["\[PartialD]", i], exp}] // ToBoxes // InterpretationBox[#,obj] &


Dx [index_,  any_] := Deriv["x"[UI[index, "lorentz"]], any]             (* Dx[\[Mu]]   = \[PartialD]^\[Mu];              \!\(\(Dx[\(,\)\(\[Mu]\)]\)\ \ \  = \ 
\*SubscriptBox[\(\[PartialD]\), \(\[Mu]\)]\) *)
Dx [,index_, any_] := Deriv["x"[LI[index, "lorentz"]], any]             (* Dx[\[Mu]]   = \[PartialD]^\[Mu];              \!\(\(Dx[\(,\)\(\[Mu]\)]\)\ \ \  = \ 
\*SubscriptBox[\(\[PartialD]\), \(\[Mu]\)]\) *)

Deriv[var_, TDot[exp__]] := TDot[DT[var], exp]
Deriv[var_, t:(_NT|_GT)] := TDot[DT[var], t]
Deriv[var_, Plus[a_, b__]] := Deriv[var, a] + Deriv[var, Plus[b]]
Deriv[var_, Times[a:(_String|_?NumericQ), b__]] := a Deriv[var, Times[b]]
Deriv[var_, (_String|_?NumericQ)] := 0

TDot[x1___, d:DT[var_], a:_NT, x2___] := TDot[x1, TakeDerivative[var, a], x2] + TDot[x1, a, d, x2]
TDot[x1___, d:DT[var_], a:_GT, x2___] := TDot[x1, TakeDerivative[var, a], x2] + If[IsGrassmann[var], -1, 1] * TDot[x1, a, d, x2]
TDot[x1___, DT[var_]] := 0

TakeDerivative["x"[UI[\[Mu]_, "lorentz"]], (t:NT|GT)[f:Field[__]|OverBar[Field[__]], any___]] := t[Dx[f], UI[\[Mu], "lorentz"], any]
TakeDerivative["x"[LI[\[Mu]_, "lorentz"]], (t:NT|GT)[f:Field[__]|OverBar[Field[__]], any___]] := t[Dx[f], LI[\[Mu], "lorentz"], any]
TakeDerivative["x"[_UI|_LI], (NT|GT)[f:Except[Field[__]|OverBar[Field[__]]], ___]] := 0
Dx[Field[n_]] := Field[Dx[n]]
Dx[OverBar[Field[n_]]] := OverBar[Field[Dx[n]]]
NT[Field[Dx[Dx[n_]]], OrderlessPatternSequence[UI[a_, "lorentz"], LI[a_, "lorentz"]], any___] := NT[Field[Dx2[n]], any]
GT[Field[Dx[Dx[n_]]], OrderlessPatternSequence[UI[a_, "lorentz"], LI[a_, "lorentz"]], any___] := GT[Field[Dx2[n]], any]

NT[f:Field[Dx[Dx[_]]], a:(UI|LI)[i_, "lorentz"], b:(UI|LI)[j_, "lorentz"], any___] /; Not[OrderedQ[{i, j}]] := NT[f, b, a, any]
GT[f:Field[Dx[Dx[_]]], a:(UI|LI)[i_, "lorentz"], b:(UI|LI)[j_, "lorentz"], any___] /; Not[OrderedQ[{i, j}]] := GT[f, b, a, any]

NT[f:Field["g"], a:(UI|LI)[i_, "lorentz"], b:(UI|LI)[j_, "lorentz"], any___] /; Not[OrderedQ[{i, j}]] := NT[f, b, a, any]


SumToList[exp_]  := If[Head[exp] === Plus, List@@exp, {exp}];

TDotOrder$Fields = <|"g"->1|>;
TDotOrder$Indices = <| |>;
TDotOrder$Params = <| |>;
TDotOrder[exp_] := Switch[exp,
  _?NumericQ, {-100, 0, exp},
  NT[Field[Dx[_]], ___], Module[{tmp=TDotOrder[exp/.Dx[a_]:>a]}, {tmp[[1]], tmp[[2]]-1, exp}],
  NT[Field["g"], ___], {5, 0, exp},
  NT[Field["\[CapitalOmega]"], ___], {3, 0, exp},
  _,    {0, 0, exp}]
SortTDot[TDot[exp__]] := TDot @@ SortBy[{exp}, TDotOrder]
SortTDot[exp_] := exp /. {f_TDot :> SortTDot[f]}

RewriteIndicesOrder[exp_] := Module[{f = TDotOrder[exp]}, {-f[[1]], f[[2;;]]}]
RewriteIndices[exp_, format_, indices_:Null] :=
  Module[{rep, old, nottoreplace, toreplace, new, rule},
    rep = exp //. TDot[a__]:>TDotList@@SortBy[{a}, RewriteIndicesOrder];
    old = Cases[rep,  format :> index, All] // Select[#, FreeQ[#, _Integer]&]&;
    nottoreplace = Counts[old] /. Association -> List // Cases[#, (p_ -> 1) :> p] &;
    toreplace = Counts[old] /. Association -> List // Cases[#, (p_ -> 2) :> p] &;
    new = If[indices===Null, If[Head[#]===OverDot, OverDot[Unique[]], Unique[]]& /@ toreplace, Select[indices, FreeQ[nottoreplace, #] &]];
    If[Length[toreplace] > Length[new], Return[exp]];
    rule = (#[[1]] -> #[[2]]) & /@ Transpose[{toreplace, new[[1 ;; Length[toreplace]]]}];
    exp /. rule];
RewriteIndices[exp_] := Module[{result},
  result = SumToList[Expand[exp]];
  result = RewriteIndices[#, (UI|LI)[index_, "lorentz"], {"\[Mu]", "\[Nu]", "\[Rho]", "\[CapitalSigma]", "\[Alpha]", "\[Beta]", "\[Gamma]"}] & /@ result;
  Total[result]
];


LowerG[exp_NT] := Which[
   MatchQ[exp[[1]], Field[Dx["g"]]], exp/.{NT[Field[Dx["g"]], x_, UI[m_,"lorentz"],UI[n_,"lorentz"]]:>Module[{a=Unique[], b=Unique[]}, -gU[m,a] Deriv["x"[x], gL[a,b]]gU[b,n]]},
   MatchQ[exp[[1]], Field[Dx[Dx[_]]]], exp/.{NT[Field[Dx[Dx[a_]]], x_, any_]:> Deriv["x"[x], LowerG[NT[Field[Dx[a]], any]]]},
   True, exp]
LowerG[exp_] := exp /. {a_NT:>LowerG[a]}


\[CapitalGamma][\[Alpha]_,\[Nu]_,\[Mu]_] := Module[{\[Beta]=Unique[]},
   TDot[(1/2)gU[\[Alpha],\[Beta]](Dx[,\[Mu],gL[\[Beta],\[Nu]]]+Dx[,\[Nu],gL[\[Beta],\[Mu]]]-Dx[,\[Beta],gL[\[Mu],\[Nu]]])]]

R[m_, n_, r_, s_] := Module[{a=Unique[]},
  Dx[,n,\[CapitalGamma][s,m,r]]-Dx[,m,\[CapitalGamma][s,n,r]]+\[CapitalGamma][a,m,r]\[CapitalGamma][s,a,n]-\[CapitalGamma][a,n,r]\[CapitalGamma][s,a,m]]
R[a_,c_] := Module[{b=Unique[]}, R[a,b,c,b]];
R[] := Module[{a=Unique[], b=Unique[]}, TDot[gU[a,b],R[a,b]]]


\[CapitalOmega] = NT[Field["\[CapitalOmega]"]];
\[CapitalOmega]inv = NT[Field["\[Omega]"]];
TDot[a___, \[CapitalOmega], b___, \[CapitalOmega]inv, c___] := TDot[a,b,c]
TDot[a___, \[CapitalOmega]inv, b___, \[CapitalOmega], c___] := TDot[a,b,c]

\[CapitalGamma]2[\[Alpha]_,\[Nu]_,\[Mu]_] := Module[{\[Beta]=Unique[]},
   TDot[(1/2)\[CapitalOmega]inv \[CapitalOmega]inv gU[\[Alpha],\[Beta]],(Dx[,\[Mu],\[CapitalOmega] \[CapitalOmega] gL[\[Beta],\[Nu]]]+Dx[,\[Nu],\[CapitalOmega] \[CapitalOmega] gL[\[Beta],\[Mu]]]-Dx[,\[Beta],\[CapitalOmega] \[CapitalOmega] gL[\[Mu],\[Nu]]])]]

R2[m_, n_, r_, s_] := Module[{a=Unique[]},
  Dx[,n,\[CapitalGamma]2[s,m,r]]-Dx[,m,\[CapitalGamma]2[s,n,r]]+\[CapitalGamma]2[a,m,r]\[CapitalGamma]2[s,a,n]-\[CapitalGamma]2[a,n,r]\[CapitalGamma]2[s,a,m]]
R2[a_,c_] := Module[{b=Unique[]}, R2[a,b,c,b]];
R2[] := Module[{a=Unique[], b=Unique[]}, TDot[\[CapitalOmega]inv, \[CapitalOmega]inv, gU[a,b],R2[a,b]]]


TDot[\[CapitalOmega],\[CapitalGamma]2[a,b,c]-\[CapitalGamma][a,b,c]]//ExpandAll// RewriteIndices


\[CapitalOmega] \[CapitalOmega] R2[] - R[] // LowerG//ExpandAll//RewriteIndices;
%//.{  NT[Field[Dx["\[Omega]"]], a_] :> -TDot[\[CapitalOmega]inv,\[CapitalOmega]inv,NT[Field[Dx["\[CapitalOmega]"]], a]]};
%//ExpandAll//RewriteIndices//SortTDot


-6\[CapitalOmega]inv gU[\[Mu],\[Nu]](Dx[,\[Mu],Dx[,\[Nu],\[CapitalOmega]]] - (Dx[,\[Rho],\[CapitalOmega]]\[CapitalGamma][\[Rho],\[Nu],\[Mu]]))//LowerG// RewriteIndices  // SortTDot//ExpandAll
