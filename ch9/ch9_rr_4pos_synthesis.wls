#!/usr/bin/env wolframscript

(* Chapter 9 spherical RR dyad synthesis for 4 specified orientations.
   Uses Eq. (9.25), bilinear forms (9.28), (9.30), center-axis cone (9.45), (9.46),
   and circling-axis cone (9.53), (9.54). *)

ClearAll["Global`*"];
$HistoryLength = 0;

prec = 70;
solvePrec = 50;
epsAxis = SetPrecision[10^-10, prec];
epsRoot = SetPrecision[10^-14, prec];
epsImag = SetPrecision[10^-18, prec];
deg = Degree;

varsG = {x, y, z};
varsW = {lam, mu, nu};

normalizeVec[v_List] := Module[{n = Sqrt[v.v]}, If[n == 0, v, v/n]];

canonicalAxis[v_List] := Module[{u, s},
  u = normalizeVec[Chop[v, epsRoot]];
  s = Which[
    u[[1]] < 0, -1,
    u[[1]] == 0 && u[[2]] < 0, -1,
    u[[1]] == 0 && u[[2]] == 0 && u[[3]] < 0, -1,
    True, 1
  ];
  Chop[s u, epsRoot]
];

dedupeScalars[list_List, tol_] := Module[{out = {}, sorted = Sort[list]},
  Do[
    If[out === {} || Abs[v - Last[out]] > tol, AppendTo[out, v]],
    {v, sorted}
  ];
  out
];

dedupeAxes[list_List, tol_] := Module[{out = {}},
  Do[
    If[AllTrue[out, Norm[# - v] > tol &], AppendTo[out, v]],
    {v, list}
  ];
  SortBy[out, {#[[1]] &, #[[2]] &, #[[3]] &}]
];

skew[v_List] := {
  {0, -v[[3]], v[[2]]},
  {v[[3]], 0, -v[[1]]},
  {-v[[2]], v[[1]], 0}
};

rotAxisAngle[axis_List, ang_] := Module[{s, k},
  s = normalizeVec[axis];
  k = skew[s];
  Chop[
    IdentityMatrix[3] + Sin[ang] k + (1 - Cos[ang]) (k.k),
    epsRoot
  ]
];

axisAngleFromRotation[a_] := Module[{phi, s},
  phi = ArcCos[Clip[(Tr[a] - 1)/2, {-1, 1}]];
  s = If[
    Abs[Sin[phi]] < 10^-22,
    {0, 0, 0},
    {
      a[[3, 2]] - a[[2, 3]],
      a[[1, 3]] - a[[3, 1]],
      a[[2, 1]] - a[[1, 2]]
    }/(2 Sin[phi])
  ];
  <|"phi" -> N[phi, 40], "S" -> N[normalizeVec[s], 40]|>
];

realRoots1D[poly_, var_Symbol] := Module[{r, rr},
  r = var /. NSolve[poly == 0, var, WorkingPrecision -> solvePrec];
  rr = Select[r, Abs[Im[#]] < epsImag &];
  dedupeScalars[Re[rr], epsRoot]
];

(* Deterministic 4-orientation task set generated from a known spherical RR motion. *)
gTrue = N[normalizeVec[{45, -32, 834}], prec];
wTrue = N[normalizeVec[{-61, 73, 31}], prec];
beta = N[{0, 35, -25, 60} deg, prec];
alpha = N[{0, -20, 30, -40} deg, prec];

A = Table[
  rotAxisAngle[gTrue, beta[[i]]] . rotAxisAngle[wTrue, alpha[[i]]],
  {i, 1, 4}
];
Arel = Table[A[[i]] . Transpose[A[[1]]], {i, 1, 4}];

exprs = Table[
  Expand[varsG . ((Arel[[i]] - IdentityMatrix[3]).varsW)],
  {i, 2, 4}
];

M = Table[
  Coefficient[exprs[[i]], #] & /@ varsW,
  {i, 1, 3}
];
M0 = Table[
  Coefficient[exprs[[i]], #] & /@ varsG,
  {i, 1, 3}
];

detCenter = Expand[Det[M]];
detCircling = Expand[Det[M0]];

(* Deterministic sampling of real center-axis directions on the cone. *)
xGrid = SetPrecision[Subdivide[-2, 2, 50], prec];
yGrid = SetPrecision[Subdivide[-2, 2, 50], prec];

chart1 = Reap[
  Do[
    rootsY = realRoots1D[Expand[detCenter /. {x -> xv, z -> 1}], y];
    Do[Sow[canonicalAxis[{xv, yv, 1}]], {yv, rootsY}],
    {xv, xGrid}
  ]
][[2]];
chart1Axes = If[chart1 === {}, {}, chart1[[1]]];

chart2 = Reap[
  Do[
    rootsZ = realRoots1D[Expand[detCenter /. {x -> xv, y -> 1}], z];
    Do[Sow[canonicalAxis[{xv, 1, zv}]], {zv, rootsZ}],
    {xv, xGrid}
  ]
][[2]];
chart2Axes = If[chart2 === {}, {}, chart2[[1]]];

chart3 = Reap[
  Do[
    rootsZ = realRoots1D[Expand[detCenter /. {y -> yv, x -> 1}], z];
    Do[Sow[canonicalAxis[{1, yv, zv}]], {zv, rootsZ}],
    {yv, yGrid}
  ]
][[2]];
chart3Axes = If[chart3 === {}, {}, chart3[[1]]];

allCenterAxes = dedupeAxes[
  Join[chart1Axes, chart2Axes, chart3Axes],
  epsAxis
];

If[allCenterAxes === {},
  Print["No real center-axis samples found on the cone."];
  Exit[1];
];

sampledCenterAxes = Take[allCenterAxes, UpTo[30]];

solveWForG[g_List] := Module[{mg, ns},
  mg = N[M /. Thread[varsG -> g], prec];
  ns = NullSpace[mg, Tolerance -> 10^-22];
  If[ns === {}, Missing["NoNullspace"], canonicalAxis[First[ns]]]
];

dyadSamples = Reap[
  Do[
    w = solveWForG[g];
    If[w =!= Missing["NoNullspace"],
      res = Table[
        N[g . ((Arel[[i]] - IdentityMatrix[3]).w), 40],
        {i, 2, 4}
      ];
      maxRes = Max[Abs[res]];
      centerRes = N[detCenter /. Thread[varsG -> g], 40];
      circlingRes = N[detCircling /. Thread[varsW -> w], 40];
      Sow[<|
        "G" -> N[g, 20],
        "W1" -> N[w, 20],
        "residuals_925" -> N[res, 30],
        "maxResidual" -> N[maxRes, 30],
        "centerConeResidual" -> N[Abs[centerRes], 30],
        "circlingConeResidual" -> N[Abs[circlingRes], 30]
      |>]
    ],
    {g, sampledCenterAxes}
  ]
][[2, 1]];

validDyads = Select[dyadSamples, #["maxResidual"] < 10^-16 &];

If[validDyads === {},
  Print["No valid dyad samples passed the residual threshold."];
  Exit[1];
];

sortedDyads = SortBy[
  validDyads,
  {
    #["G"][[1]] &,
    #["G"][[2]] &,
    #["G"][[3]] &
  }
];

monomialsCenter = {
  y^3, x y^2, z y^2, x^2 y, x y z, z^2 y, x^3, x^2 z, x z^2, z^3
};
centerCoeffs = AssociationThread[
  ToString /@ monomialsCenter,
  N[Coefficient[detCenter, #], 30] & /@ monomialsCenter
];

monomialsCircling = {
  mu^3, lam mu^2, nu mu^2, lam^2 mu, lam nu mu, nu^2 mu,
  lam^3, lam^2 nu, lam nu^2, nu^3
};
circlingCoeffs = AssociationThread[
  ToString /@ monomialsCircling,
  N[Coefficient[detCircling, #], 30] & /@ monomialsCircling
];

axisAngleTable = Association@Table[
  "A1" <> ToString[i] -> axisAngleFromRotation[Arel[[i]]],
  {i, 2, 4}
];

resultAssoc = <|
  "precisionDigits" -> prec,
  "solvePrecisionDigits" -> solvePrec,
  "taskData" -> <|
    "generationDescription" -> "Synthetic 4-orientation task generated from known RR dyad (deterministic).",
    "G_true" -> N[gTrue, 20],
    "W1_true" -> N[wTrue, 20],
    "beta_deg" -> {0, 35, -25, 60},
    "alpha_deg" -> {0, -20, 30, -40},
    "orientations_Ai" -> Association@Table[
      "M" <> ToString[i] -> N[A[[i]], 20],
      {i, 1, 4}
    ],
    "relative_A1i" -> Association@Table[
      "A1" <> ToString[i] -> N[Arel[[i]], 20],
      {i, 2, 4}
    ],
    "relativeAxisAngle" -> axisAngleTable
  |>,
  "centerCone" -> <|
    "equation" -> ToString[InputForm[detCenter]],
    "coefficients_946_style" -> centerCoeffs
  |>,
  "circlingCone" -> <|
    "equation" -> ToString[InputForm[detCircling]],
    "coefficients_954_style" -> circlingCoeffs
  |>,
  "sampling" -> <|
    "numSampledCenterAxes" -> Length[sampledCenterAxes],
    "numValidDyads" -> Length[sortedDyads]
  |>,
  "dyadSamples" -> MapIndexed[
    Join[<|"sampleId" -> First[#2]|>, #1] &,
    sortedDyads
  ]
|>;

csvHeader = {
  "sampleId",
  "Gx", "Gy", "Gz",
  "W1x", "W1y", "W1z",
  "res2", "res3", "res4",
  "maxResidual",
  "centerConeResidual",
  "circlingConeResidual"
};

csvRows = Prepend[
  Map[
    {
      #["sampleId"],
      #["G"][[1]], #["G"][[2]], #["G"][[3]],
      #["W1"][[1]], #["W1"][[2]], #["W1"][[3]],
      #["residuals_925"][[1]], #["residuals_925"][[2]], #["residuals_925"][[3]],
      #["maxResidual"],
      #["centerConeResidual"],
      #["circlingConeResidual"]
    } &,
    resultAssoc["dyadSamples"]
  ],
  csvHeader
];

Export["ch9_rr_4pos_results.json", resultAssoc, "JSON"];
Export["ch9_rr_4pos_samples.csv", csvRows];
Export["ch9_rr_4pos_center_cone.txt", ToString[InputForm[detCenter]], "Text"];
Export["ch9_rr_4pos_circling_cone.txt", ToString[InputForm[detCircling]], "Text"];

Print["Generated files:"];
Print[" - ch9_rr_4pos_results.json"];
Print[" - ch9_rr_4pos_samples.csv"];
Print[" - ch9_rr_4pos_center_cone.txt"];
Print[" - ch9_rr_4pos_circling_cone.txt"];
Print["Valid dyad samples: ", Length[resultAssoc["dyadSamples"]]];
