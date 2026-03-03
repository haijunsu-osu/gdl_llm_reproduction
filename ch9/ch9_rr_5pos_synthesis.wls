#!/usr/bin/env wolframscript

(* Chapter 9 spherical RR dyad synthesis for 5 specified orientations.
   Implements Eq. (9.25), (9.55), (9.56), (9.57), (9.58), (9.59). *)

ClearAll["Global`*"];
$HistoryLength = 0;

prec = 80;
solvePrec = 60;
epsImag = SetPrecision[10^-16, prec];
epsRoot = SetPrecision[10^-14, prec];
epsAxis = SetPrecision[10^-10, prec];
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

sortComplex[list_List] := SortBy[list, {Re[#] &, Im[#] &}];

dedupeByDistance[list_List, tol_] := Module[{out = {}},
  Do[
    If[AllTrue[out, Abs[# - v] > tol &], AppendTo[out, v]],
    {v, list}
  ];
  out
];

groupRootsWithMultiplicity[roots_List, tol_] := Module[
  {groups = {}, added = False, rep, i},
  Do[
    added = False;
    For[i = 1, i <= Length[groups], i++,
      rep = groups[[i, "representative"]];
      If[Abs[r - rep] <= tol,
        groups[[i, "roots"]] = Append[groups[[i, "roots"]], r];
        groups[[i, "representative"]] = Mean[groups[[i, "roots"]]];
        added = True;
        Break[];
      ];
    ];
    If[!added,
      groups = Append[groups, <|"representative" -> r, "roots" -> {r}|>];
    ],
    {r, roots}
  ];
  MapIndexed[
    <|
      "groupId" -> First[#2],
      "representative" -> #["representative"],
      "multiplicity" -> Length[#["roots"]],
      "roots" -> #["roots"],
      "isReal" -> (Abs[Im[#["representative"]]] < epsImag)
    |>&,
    groups
  ]
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
    Abs[Sin[phi]] < 10^-24,
    {0, 0, 0},
    {
      a[[3, 2]] - a[[2, 3]],
      a[[1, 3]] - a[[3, 1]],
      a[[2, 1]] - a[[1, 2]]
    }/(2 Sin[phi])
  ];
  <|"phi" -> N[phi, 40], "S" -> N[normalizeVec[s], 40]|>
];

(* Deterministic synthetic 5-orientation task generated from a known spherical RR dyad. *)
gTrue = N[normalizeVec[{51, -37, 910}], prec];
wTrue = N[normalizeVec[{-73, 64, 29}], prec];
beta = N[{0, 28, -33, 47, -59} deg, prec];
alpha = N[{0, -18, 24, -31, 42} deg, prec];

A = Table[
  rotAxisAngle[gTrue, beta[[i]]] . rotAxisAngle[wTrue, alpha[[i]]],
  {i, 1, 5}
];
Arel = Table[A[[i]] . Transpose[A[[1]]], {i, 1, 5}];
i3 = IdentityMatrix[3];

(* Eq. (9.25): G . (A1i - I) . W1 = 0, i=2..5 *)
eq925 = Table[
  Expand[varsG . ((Arel[[i]] - i3).varsW)],
  {i, 2, 5}
];

(* Eq. (9.55): 4x3 matrix M from coefficients in lambda, mu, nu *)
M = Table[
  Coefficient[eq925[[r]], #] & /@ varsW,
  {r, 1, 4}
];

(* Eq. (9.56): minors R_j by removing row 5-j *)
R = Table[
  Expand[Det[Delete[M, 5 - j]]],
  {j, 1, 4}
];

(* Eq. (9.57): dehomogenize z=1 and collect powers of y *)
Rxy = Expand[# /. z -> 1] & /@ R;

dRow[poly_] := Coefficient[poly, y, #] & /@ {3, 2, 1, 0};  (* {d_j0,d_j1,d_j2,d_j3} *)
Dmat = dRow /@ Rxy;

(* Eq. (9.59): P(x) = det(D) *)
P = Expand[Det[Dmat]];
pDegree = Exponent[P, x];

If[pDegree > 6 || pDegree < 1,
  Print["Unexpected elimination degree for P(x): ", pDegree];
];

pCoeffs = Association@Table[
  "x^" <> ToString[k] -> N[Coefficient[P, x, k], 40],
  {k, pDegree, 0, -1}
];

rootsAllRaw = x /. NSolve[P == 0, x, WorkingPrecision -> solvePrec];
rootsAll = sortComplex[rootsAllRaw];

rootGroups = groupRootsWithMultiplicity[rootsAll, 10^-12];
realRootGroups = Select[rootGroups, #["isReal"] &];

(* For each real root x_i, solve for y_i from (9.58)/(9.57), then recover G_i and W1_i. *)
candidateSolutions = Reap[
  Do[
    xi = Re[group["representative"]];

    Dxi = N[Dmat /. x -> xi, prec];
    nsD = NullSpace[Dxi, Tolerance -> 10^-24];
    yFromNull = If[
      nsD === {} || Abs[First[nsD][[4]]] < 10^-24,
      Missing["NoNullY"],
      Re[First[nsD][[3]]/First[nsD][[4]]]
    ];

    yRootsRaw = y /. NSolve[(Rxy[[1]] /. x -> xi) == 0, y, WorkingPrecision -> solvePrec];
    yRootsReal = Select[yRootsRaw, Abs[Im[#]] < epsImag &];
    yCandidates = dedupeByDistance[Re[yRootsReal], 10^-10];
    If[yFromNull =!= Missing["NoNullY"],
      yCandidates = dedupeByDistance[Append[yCandidates, yFromNull], 10^-10];
    ];

    Do[
      rVals = N[(Rxy /. {x -> xi, y -> yi}), 40];
      maxR = Max[Abs[rVals]];
      If[maxR < 10^-9,
        gHom = {xi, yi, 1};
        gAxis = canonicalAxis[gHom];

        (* As in Section 9.4.6.1, choose two constraints (9.25) to solve for lambda,mu,nu. *)
        m2 = N[(M /. Thread[varsG -> gHom])[[1 ;; 2, All]], prec];
        nsW = NullSpace[m2, Tolerance -> 10^-24];
        If[nsW =!= {},
          wAxis = canonicalAxis[First[nsW]];
          res925 = Table[
            N[gAxis . ((Arel[[i]] - i3).wAxis), 40],
            {i, 2, 5}
          ];
          maxRes925 = Max[Abs[res925]];
          fullM = N[M /. Thread[varsG -> gHom], prec];
          resMAll = N[fullM.wAxis, 40];
          maxResMAll = Max[Abs[resMAll]];

          Sow[<|
            "rootGroupId" -> group["groupId"],
            "rootRepresentative_x" -> N[group["representative"], 40],
            "rootMultiplicity" -> group["multiplicity"],
            "x" -> N[xi, 40],
            "y" -> N[yi, 40],
            "RxyResiduals" -> N[rVals, 30],
            "maxRxyResidual" -> N[maxR, 30],
            "G_homogeneous" -> N[gHom, 30],
            "G_axis" -> N[gAxis, 30],
            "W1_axis" -> N[wAxis, 30],
            "eq925Residuals" -> N[res925, 30],
            "maxEq925Residual" -> N[maxRes925, 30],
            "matrixResiduals_MW" -> N[resMAll, 30],
            "maxMatrixResidual_MW" -> N[maxResMAll, 30]
          |>];
        ];
      ],
      {yi, yCandidates}
    ],
    {group, realRootGroups}
  ]
][[2]];

allCandidates = If[candidateSolutions === {}, {}, candidateSolutions[[1]]];

(* Deduplicate geometric solutions by G axis and W1 axis. *)
sortedCandidates = SortBy[
  allCandidates,
  {
    #["x"] &,
    #["y"] &
  }
];

uniqueSolutions = Reap[
  seen = {};
  Do[
    key = Join[s["G_axis"], s["W1_axis"]];
    If[AllTrue[seen, Norm[# - key] > epsAxis &],
      AppendTo[seen, key];
      Sow[s];
    ],
    {s, sortedCandidates}
  ]
];
finalSolutions = If[uniqueSolutions[[2]] === {}, {}, uniqueSolutions[[2, 1]]];

validatedSolutions = Select[finalSolutions, #["maxEq925Residual"] < 10^-12 &];

axisAngleTable = Association@Table[
  "A1" <> ToString[i] -> axisAngleFromRotation[Arel[[i]]],
  {i, 2, 5}
];

rootsAssoc = MapIndexed[
  <|
    "rootId" -> First[#2],
    "realPart" -> N[Re[#1], 40],
    "imagPart" -> N[Im[#1], 40]
  |>&,
  rootsAll
];

rootGroupsAssoc = Map[
  <|
    "groupId" -> #["groupId"],
    "representativeRealPart" -> N[Re[#["representative"]], 40],
    "representativeImagPart" -> N[Im[#["representative"]], 40],
    "multiplicity" -> #["multiplicity"],
    "isReal" -> #["isReal"]
  |>&,
  rootGroups
];

resultAssoc = <|
  "precisionDigits" -> prec,
  "solvePrecisionDigits" -> solvePrec,
  "taskData" -> <|
    "generationDescription" -> "Synthetic deterministic 5-orientation task generated from a known spherical RR dyad.",
    "G_true" -> N[gTrue, 30],
    "W1_true" -> N[wTrue, 30],
    "beta_deg" -> {0, 28, -33, 47, -59},
    "alpha_deg" -> {0, -18, 24, -31, 42},
    "orientations_Ai" -> Association@Table[
      "M" <> ToString[i] -> N[A[[i]], 30],
      {i, 1, 5}
    ],
    "relative_A1i" -> Association@Table[
      "A1" <> ToString[i] -> N[Arel[[i]], 30],
      {i, 2, 5}
    ],
    "relativeAxisAngle" -> axisAngleTable
  |>,
  "elimination" -> <|
    "M_955" -> ToString[InputForm[M]],
    "Rj_956" -> Association@Table[
      "R" <> ToString[j] -> ToString[InputForm[R[[j]]]],
      {j, 1, 4}
    ],
    "Rj_dehomogenized_957" -> Association@Table[
      "R" <> ToString[j] -> ToString[InputForm[Rxy[[j]]]],
      {j, 1, 4}
    ],
    "D_958" -> ToString[InputForm[Dmat]],
    "P_959_degree" -> pDegree,
    "P_959_equation" -> ToString[InputForm[P]],
    "P_959_coefficients" -> pCoeffs
  |>,
  "roots" -> <|
    "allRootsSorted" -> rootsAssoc,
    "rootGroupsWithMultiplicity" -> rootGroupsAssoc
  |>,
  "solutions" -> <|
    "numCandidatesFromRealRoots" -> Length[allCandidates],
    "numUniqueGeometricSolutions" -> Length[finalSolutions],
    "numValidatedSolutions" -> Length[validatedSolutions],
    "validated" -> MapIndexed[
      Join[<|"solutionId" -> First[#2]|>, #1] &,
      validatedSolutions
    ]
  |>
|>;

rootsCsvHeader = {"rootId", "realPart", "imagPart"};
rootsCsvRows = Prepend[
  (List @@ # & /@ (Normal /@ rootsAssoc)),
  rootsCsvHeader
];

groupCsvHeader = {"groupId", "representativeRealPart", "representativeImagPart", "multiplicity", "isReal"};
groupCsvRows = Prepend[
  (List @@ # & /@ (Normal /@ rootGroupsAssoc)),
  groupCsvHeader
];

solHeader = {
  "solutionId", "rootGroupId", "rootMultiplicity", "x", "y",
  "Gx", "Gy", "Gz", "W1x", "W1y", "W1z",
  "res2", "res3", "res4", "res5", "maxEq925Residual"
};
solRows = Prepend[
  Map[
    {
      #["solutionId"],
      #["rootGroupId"],
      #["rootMultiplicity"],
      #["x"],
      #["y"],
      #["G_axis"][[1]], #["G_axis"][[2]], #["G_axis"][[3]],
      #["W1_axis"][[1]], #["W1_axis"][[2]], #["W1_axis"][[3]],
      #["eq925Residuals"][[1]], #["eq925Residuals"][[2]],
      #["eq925Residuals"][[3]], #["eq925Residuals"][[4]],
      #["maxEq925Residual"]
    } &,
    resultAssoc["solutions", "validated"]
  ],
  solHeader
];

Export["ch9_rr_5pos_results.json", resultAssoc, "JSON"];
Export["ch9_rr_5pos_roots.csv", rootsCsvRows];
Export["ch9_rr_5pos_root_groups.csv", groupCsvRows];
Export["ch9_rr_5pos_solutions.csv", solRows];
Export["ch9_rr_5pos_Px.txt", ToString[InputForm[P]], "Text"];

Print["Generated files:"];
Print[" - ch9_rr_5pos_results.json"];
Print[" - ch9_rr_5pos_roots.csv"];
Print[" - ch9_rr_5pos_root_groups.csv"];
Print[" - ch9_rr_5pos_solutions.csv"];
Print[" - ch9_rr_5pos_Px.txt"];
Print["P(x) degree: ", pDegree];
Print["Validated solutions: ", Length[resultAssoc["solutions", "validated"]]];
