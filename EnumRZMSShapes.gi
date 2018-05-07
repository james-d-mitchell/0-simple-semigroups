# This code requires the images package
LoadPackage("images");

###############################################################################
# Permutations and operations on sets represenating (regular) binary matrices
###############################################################################
# Calculates the perm group which acts on these sets (via OnSets) by
# transposition, row permutations and column permutations
ActingBinaryMatrixGroup := function(nr_rows, nr_cols)
  local out, temp, a, i;

  if nr_rows = 1 and nr_cols = 1 then
    return Group(PermList([1]));
  fi;

  out := [];
  if nr_rows = nr_cols then
    temp := [];
    for i in [1 .. nr_rows] do
      Append(temp, List([1 .. nr_cols], j -> i + (j - 1) * nr_cols));
    od;
    Add(out, PermList(temp));
  fi;

  if nr_rows > 1 then
    Add(out, PermList(
      Concatenation([nr_cols + 1 .. 2 * nr_cols], [1 .. nr_cols])));
    if nr_rows > 2 then
      Add(out, PermList(
        Concatenation([(nr_rows - 1) * nr_cols + 1 .. nr_rows * nr_cols],
                      [1 .. (nr_rows - 1) * nr_cols])));
    fi;
  fi;

  if nr_cols > 1 then
    temp := Concatenation([2, 1], [3 .. nr_cols]);
    Add(out, PermList(
      Concatenation(temp + nr_cols, temp,
        Concatenation(List([3 .. nr_rows], a -> temp + nr_cols * (a - 1))))));
    if nr_cols > 2 then
      temp := Concatenation([2 .. nr_cols], [1]);
      Add(out, PermList(
        Concatenation(List([1 .. nr_rows], a -> temp + nr_cols * (a - 1)))));
    fi;
  fi;

  return Group(out);
end;

###############################################################################
# This permutation sends (via OnSets) a set representing an m x n matrix to its
# transpose as a set representing a n x m matrix.
TranspositionPermutation := function(nr_rows, nr_cols)
  return PermList(Concatenation(
           List([1 .. nr_rows], i ->
             List([1 .. nr_cols], j -> (j - 1) * nr_rows) + i)));
end;

###############################################################################
# This returns the binary matrix represented by the set <set>. Note this
# depends on the dimensions!
SetToBinaryMat := function(set, nr_rows, nr_cols)
  local func, out, i, row;
  func := function(b)
    if b then
      return 1;
    fi;
    return 0;
  end;
  out := EmptyPlist(nr_rows);
  for i in [1 .. nr_rows] do
    row := List([(i - 1) * nr_cols + 1 .. i * nr_cols], a -> func(a in set));
    Add(out, ShallowCopy(row));
  od;
  return out;
end;

BinaryMatToSet := function(mat)
  local out, i, j;
  out := [];
  for i in [1 .. Length(mat)] do
    for j in [1 .. Length(mat[1])] do
      if mat[i][j] = 1 then
        AddSet(out, (i - 1) * Length(mat[1]) + j);
      fi;
    od;
  od;
  return out;
end;

BooleanMatToSet := function(mat)
  local out, i, j;
  out := [];
  for i in [1 .. Length(mat)] do
    for j in [1 .. Length(mat[1])] do
      if mat[i][j] then
        AddSet(out, (i - 1) * Length(mat[1]) + j);
      fi;
    od;
  od;
  return out;
end;

###############################################################################
# This returns a transformation which relabels the elements of the set so that
# it represents the upper left corner of a matrix with more rows and/or columns
EmbedSetTransformation := function(nr_rows, nr_cols, new_rows, new_cols)
  local trans;
  trans := Concatenation(List([1 .. nr_rows], i -> List([1 .. nr_cols], j ->
           j + (i - 1) * new_cols)));
  Append(trans, List([1 .. (new_cols - nr_cols) * (nr_rows - 1)], a -> 1));
  return TransformationList(trans);
end;

###############################################################################
# Returns the extensions for a n-1 x n-1 01matrix (as a set) to a n x n matrix.
SquareExtensions := function(n)
  local last_col, last_row, ext1, ext2, comb1, comb2, c;
  last_col := Set([1 .. n - 1], i -> i * n);
  last_row := Set([1 .. n - 1], i -> n * (n - 1) + i);

  ext1 := Combinations(Union(last_col, last_row));
  for c in ext1 do
    AddSet(c, n ^ 2);
  od;

  comb1 := Combinations(last_col);
  Remove(comb1, 1);
  comb2 := Combinations(last_row);
  Remove(comb2, 1);

  ext2 := Cartesian(comb1, comb2);
  Apply(ext2, Union);

  return Union(ext1, ext2);
end;

# Returns the extensions for a m-1 x n 01matirx (as a set) to a m x n matrix.
RowExtensions := function(nr_rows, nr_cols)
  local last_row, cmb;
  last_row := Set([1 .. nr_cols], i -> nr_cols * (nr_rows - 1) + i);

  cmb := Combinations(last_row);
  Remove(cmb, 1); # Remove the empty row.

  return cmb;
end;

###############################################################################
###############################################################################
# Finding the binary matrix shapes
###############################################################################
###############################################################################
# Solution to two row case.
BinaryMatrixShapesWithTwoRows := function(nr_cols)
  local row1, row2, out, i, j, opp;
  out := [];
  row1 := [1 .. nr_cols];
  for j in [1 .. nr_cols] do
    Add(out, Union(row1, [1 .. j] + nr_cols));
  od;

  i := nr_cols - 1;
  while 2 * i >= nr_cols do
    row1 := [1 .. i];
    opp := [1 .. nr_cols];
    SubtractSet(opp, row1);
    opp := opp + nr_cols;
    for j in [0 .. 2 * i - nr_cols] do
      row2 := Union(opp, [1 .. j] + nr_cols);
      Add(out, Union(row1, row2));
    od;
    i := i - 1;
  od;
  return out;
end;

###############################################################################
# Wrapper for solved cases.
BinaryMatrixShapesSolvedCases := function(nr_rows, nr_cols)
  if nr_rows = 1 then
    return [[1 .. nr_cols]];
  elif nr_cols = 1 then
    return List([1 .. nr_rows], a -> [a]);
  elif nr_rows = 2 then
    return BinaryMatrixShapesWithTwoRows(nr_cols);
  fi;
end;

###############################################################################
###############################################################################
# IO methods used whilst seraching as well as to save results.
###############################################################################
###############################################################################
BinaryMatrixShapesReadFile := function(nr_rows, nr_cols)
  local prefix, name, file, out;
  prefix := "/Users/crussell/Desktop/shapes/shapes-";
  name := Concatenation(prefix, String(nr_rows), "x", String(nr_cols));
  file := IO_File(name, "r");
  if file = fail then
    return fail;
  fi;
  out := IO_ReadLines(file);
  IO_Close(file);
  Apply(out, a -> EvalString(Chomp(a)));
  return out;
end;

BinaryMatrixShapesWriteFile := function(nr_rows, nr_cols, data)
  local prefix, name, file;
  prefix := "/Users/crussell/Desktop/shapes/shapes-";
  name := Concatenation(prefix, String(nr_rows), "x", String(nr_cols));
  file := IO_File(name, "w");
  IO_WriteLines(file, data);
  IO_Close(file);
end;

###############################################################################
###############################################################################
# Main method:
# This finds shapes by building up from lower dimension cases.
###############################################################################
###############################################################################
BinaryMatrixShapesByDim := function(arg)
  local nr_rows, nr_cols, RW, data, out, trans, embed, G, extensions, pos, temp,
  d, e;

  if Length(arg) = 2 then
    nr_rows := arg[1];
    nr_cols := arg[2];
    RW := false;
  elif Length(arg) = 3 then
    nr_rows := arg[1];
    nr_cols := arg[2];
    RW := arg[3];
  else
    Error("BinaryMatrixShapesByDim: invalid arguments");
  fi;

  # If we have calculated this case before then return those results
  if RW then
    data := BinaryMatrixShapesReadFile(nr_rows, nr_cols);
    if data <> fail then
      return data;
    fi;
  fi;

  # These are cases which we can describe the solution to explicitly
  if nr_rows <= 2 or nr_cols = 1 then
    return BinaryMatrixShapesSolvedCases(nr_rows, nr_cols);
  fi;

  # The m x n case is deduced from the n x m case.
  if nr_rows < nr_cols or nr_cols = 2 then
    out := BinaryMatrixShapesByDim(nr_cols, nr_rows);
    trans := TranspositionPermutation(nr_rows, nr_cols);
    Apply(out, a -> OnSets(a, trans));

  # In square matrix case, build up from previous square case
  elif nr_rows = nr_cols then
    data := BinaryMatrixShapesByDim(nr_rows - 1, nr_cols - 1);
    embed := EmbedSetTransformation(nr_rows - 1, nr_cols - 1, nr_rows, nr_cols);
    Apply(data, d -> OnSets(d, embed));
    G := ActingBinaryMatrixGroup(nr_rows, nr_cols);
    extensions := SquareExtensions(nr_rows);
    out := [];
    for d in data do
      pos := Position(data, d);
      for e in extensions do
        Print("\c Rows = ", nr_rows, "; Cols = ", nr_cols, "; Data pos = ", pos,
        "/", Size(data), "; Extensions pos = ", Position(extensions, e), "/",
        Size(extensions), "\r");
        temp := ShallowCopy(d);
        UniteSet(temp, e);
        AddSet(out, CanonicalImage(G, temp, OnSets));
      od;
    od;

  # If more rows than columns then build from a square matrix by adding rows.
  elif nr_rows > nr_cols then
    data := BinaryMatrixShapesByDim(nr_rows - 1, nr_cols);
    G := ActingBinaryMatrixGroup(nr_rows, nr_cols);
    extensions := RowExtensions(nr_rows, nr_cols);
    Remove(extensions, nr_cols); # this omits the extension of a full row of 1's
    out := EmptyPlist(Size(extensions) * Size(data));
    for d in data do
      pos := Position(data, d);
      for e in extensions do
        Print("\c Rows = ", nr_rows, "; Cols = ", nr_cols, "; Data pos = ", pos,
        "/", Size(data), "; Extensions pos = ", Position(extensions, e), "/",
        Size(extensions), "\r");
        temp := ShallowCopy(d);
        UniteSet(temp, e);
        AddSet(out, CanonicalImage(G, temp, OnSets));
      od;
    od;
    AddSet(out, [1 .. nr_rows * nr_cols]); # this is the only case we miss with
                                           # earlier omitted extension
  fi;

  if RW then
    BinaryMatrixShapesWriteFile(nr_rows, nr_cols, out);
  fi;
  return out;
end;

###############################################################################

CanonicalizeBinaryMatrixSets := function(nr_rows, nr_cols, sets)
  local G, out, s;
  G := ActingBinaryMatrixGroup(nr_rows, nr_cols);
  out := [];
  for s in sets do
      AddSet(out, CanonicalImage(G, s, OnSets));
  od;
  return out;
end;

BinaryMatrixByDimTest := function(nr_rows, nr_cols, extensions)
  local data, G, out, pos, temp, d, e;
  data := BinaryMatrixShapesByDim(nr_rows - 1, nr_cols);
  G := ActingBinaryMatrixGroup(nr_rows, nr_cols);
  out := EmptyPlist(Size(extensions) * Size(data));
  for d in data do
    pos := Position(data, d);
    for e in extensions do
      Print("\c Rows = ", nr_rows, "; Cols = ", nr_cols, "; Data pos = ", pos,
      "/", Size(data), "; Extensions pos = ", Position(extensions, e), "/",
      Size(extensions), "\r");
      temp := ShallowCopy(d);
      UniteSet(temp, e);
      AddSet(out, CanonicalImage(G, temp, OnSets));
    od;
  od;
  return out;
end;

###############################################################################
# Filtering extensions
###############################################################################
# RowArrays represent a binary matrix by a subset of
# [1 .. nr_rows * (2 ^ nr_cols - 1)] where each integer represents a row in the
# row space of binary rows with nr_cols columns as well as a number in [1 ..
# nr_rows] defining how often that row occurs in the binary matrix.
# This is appropriate since we don't care about the order of the rows in the
# matrix - we want them to always be in a certain order!
_RowArrayActingGroup := function(nr_rows, nr_cols, row_space)
  local size_rs, g1, g2, g3, g4;

  size_rs := 2 ^ nr_cols - 1;
  g1 := PermList(List([1 .. 2 ^ nr_cols - 1],
        i -> Position(row_space, Permuted(row_space[i], (1, 2)))));
  g2 := PermList(List([1 .. 2 ^ nr_cols - 1],
        i -> Position(row_space, Permuted(row_space[i],
                      PermList(Concatenation([2 .. nr_cols], [1]))))));

  g3 := PermList(List([1 .. nr_rows * size_rs],
        a -> (((a + nr_rows - (1 + (a - 1) mod nr_rows)) / nr_rows) ^ g1) *
             nr_rows - nr_rows + 1 + (a - 1) mod nr_rows));
  g4 := PermList(List([1 .. nr_rows * size_rs],
        a -> (((a + nr_rows - (1 + (a - 1) mod nr_rows)) / nr_rows) ^ g1) *
             nr_rows - nr_rows + 1 + (a - 1) mod nr_rows));

  return Group(g3, g4);
end;

# Turn m x n representations into all possible (m + 1) x n extensions by
# non-zero rows
_RowArrayExtensions := function(array, nr_rows, nr_cols)
  local out, nr, rows, m, tmp, tmp2, x, i;
  out := [];
  nr := [];
  rows := [];
  for x in array do
    m := 1 + (x - 1) mod nr_rows;
    Add(nr, m);
    Add(rows, (x + nr_rows - m) / nr_rows);
  od;
  tmp := List([1 .. Size(rows)], i -> (rows[i] - 1) * (nr_rows + 1) + nr[i]);
  for i in [1 .. 2 ^ nr_cols - 1] do
    if not i in rows then
      tmp2 := ShallowCopy(tmp);
      Add(tmp2, 1 + (i - 1) * (nr_rows + 1));
      Add(out, tmp2);
    fi;
  od;
  for i in [1 .. Size(rows)] do
    tmp2 := ShallowCopy(tmp);
    tmp2[i] := tmp2[i] + 1;
    Add(out, tmp2);
  od;
  Perform(out, Sort);
  return out;
end;

_SetToRowArray := function(set, nr_rows, nr_cols, row_space)
  local size_rs, out, counts, i;
  size_rs := 2 ^ nr_cols - 1;
  out := List(SetToBinaryMat(set, nr_rows, nr_cols),
              a -> Position(row_space, a));
  counts := EmptyPlist(size_rs);
  for i in out do
    if IsBound(counts[i]) then
      counts[i] := counts[i] + 1;
    else
      counts[i] := 1;
    fi;
  od;
  out := [];
  for i in [1 .. size_rs] do
    if IsBound(counts[i]) then
      Add(out, (i - 1) * nr_rows + counts[i]);
    fi;
  od;
  return out;
end;

_RowArrayToSet := function(array, nr_rows, nr_cols, row_space)
  local ret, out, count, value, row, i, k;
    ret := List(array, a -> 1 + (a - 1) mod nr_rows);
    out := EmptyPlist(nr_rows);
    count := 0;
    for i in [1 .. Size(array)] do
      value := 1 + (array[i] - ret[i]) / (nr_rows);
      row := row_space[value];
      for k in [1 .. ret[i]] do
        out[nr_rows - count] := row;
        count := count + 1;
      od;
    od;
  return BinaryMatToSet(out);
end;

# Attempt to speed up search by first filtering extensions. To do this we
# canonicalize w.r.t rows (and assuming columns then go into dictionary order).
# data is the data from the nr_rows - 1 x nr_cols case.
_RowArrayByDim := function(nr_rows, nr_cols, data, row_space)
  local G, out, count, extensions, d, e;
  if nr_rows <> nr_cols then
    G := _RowArrayActingGroup(nr_rows, nr_cols, row_space);
    out := [];
    count := 1;
    for d in data do
      extensions := _RowArrayExtensions(d, nr_rows - 1, nr_cols);
      for e in extensions do
        Print("\c Rows = ", nr_rows, "; Cols = ", nr_cols,
        "; Row Arrays Data pos = ", count, "/", Size(data), "\r");
        AddSet(out, CanonicalImage(G, e, OnSets));
      od;
      count := count + 1;
    od;
  else
    Error("_RowArrayByDim: there should be unequal numbers of rows and ",
          "columns,");
  fi;
  return out;
end;

BinaryMatrixShapesByDimFaster := function(nr_rows, nr_cols)
  local row_space, data, out;
    row_space := Cartesian(List([1 .. nr_cols], a -> [0, 1]));
    Remove(row_space, 1);
    data := BinaryMatrixShapesByDim(nr_rows - 1, nr_cols);
    Apply(data, d -> _SetToRowArray(d, nr_rows - 1, nr_cols, row_space));
    out := _RowArrayByDim(nr_rows, nr_cols, data, row_space);
    Apply(out, d -> _RowArrayToSet(d, nr_rows, nr_cols, row_space));
    return CanonicalizeBinaryMatrixSets(nr_rows, nr_cols, out);
end;

###############################################################################
###############################################################################
# Viewing matrices:
###############################################################################
###############################################################################
PrintMat := function(m)
  local r;
  for r in m do
    Print(r, "\n");
  od;
  Print("\n");
end;

###############################################################################
# Notes:
###############################################################################
# -Use CanonicalImage with OnSets. [[1,0],[1,1]] = {1,3,4} etc.
# -Build row by row.
# -Use tricks to reduce search space, such as rows need to have less 1's than
# above row. Only do this where worthwhile!
# -Could use BlissCanonicalDigraph with m * n two coloured clique graph

###############################################################################
# Lovelace implementation:
###############################################################################
_ParallelBinaryMatrixShapesReadFile := function(nr_rows, nr_cols)
  local prefix, name, file, out;
  prefix := "/circa/home/cr66/shapes/shapes-";
  name := Concatenation(prefix, String(nr_rows), "x", String(nr_cols));
  file := IO_File(name, "r");
  if file = fail then
    return fail;
  fi;
  out := IO_ReadLines(file);
  IO_Close(file);
  Apply(out, a -> EvalString(Chomp(a)));
  return out;
end;

_ParallelBinaryMatrixShapesWriteFile := function(nr_rows, nr_cols, data,
                                                 fork_nr)
  local prefix, name, file;
  prefix := "/circa/home/cr66/shapes/shapes-";
  name := Concatenation(prefix, String(nr_rows), "x", String(nr_cols), "-fork",
          String(fork_nr));
  file := IO_File(name, "w");
  IO_WriteLines(file, data);
  IO_Close(file);
end;

# This doesn't work for nr_rows <> nr_cols!!!!!!!!
_ParallelBinaryMatrixShapesByDim := function(nr_rows, nr_cols, fork_nr,
                                             nr_forks)
  local data, range, embed, G, extensions, out, temp, i, e;

  # If we have calculated this case before then return those results
  data := _ParallelBinaryMatrixShapesReadFile(nr_rows, nr_cols);
  if data <> fail then
    return data;
  fi;

  # In square matrix case, build up from previous square case
  if nr_rows = nr_cols then
    data := _ParallelBinaryMatrixShapesReadFile(nr_rows - 1, nr_cols - 1);
    range := Int(Floor(Float(Size(data) / nr_forks)));
    if fork_nr = nr_forks then
      range := [1 + range * (nr_forks - 1) .. Size(data)];
    else
      range := [1 + range * (fork_nr - 1) .. range * fork_nr];
    fi;

    embed := EmbedSetTransformation(nr_rows - 1, nr_cols - 1, nr_rows, nr_cols);
    Apply(data, d -> OnSets(d, embed));
    G := ActingBinaryMatrixGroup(nr_rows, nr_cols);
    extensions := SquareExtensions(nr_rows);
    out := [];
    for i in range do
      for e in extensions do
        temp := ShallowCopy(data[i]);
        UniteSet(temp, e);
        AddSet(out, CanonicalImage(G, temp, OnSets));
      od;
    od;

  # If more rows than columns then build from a square matrix by adding rows.
  elif nr_rows > nr_cols then
    data := _ParallelBinaryMatrixShapesReadFile(nr_rows - 1, nr_cols);
    range := Int(Floor(Float(Size(data) / nr_forks)));
    if fork_nr = nr_forks then
      range := [1 + range * (nr_forks - 1) .. Size(data)];
    else
      range := [1 + range * (fork_nr - 1) .. range * fork_nr];
    fi;

    G := ActingBinaryMatrixGroup(nr_rows, nr_cols);
    extensions := RowExtensions(nr_rows, nr_cols);
    out := EmptyPlist(Size(extensions) * Size(data));
    for i in range do
      for e in extensions do
        temp := ShallowCopy(data[i]);
        UniteSet(temp, e);
        AddSet(out, CanonicalImage(G, temp, OnSets));
      od;
    od;
  else
    return fail;
  fi;

  _ParallelBinaryMatrixShapesWriteFile(nr_rows, nr_cols, out, fork_nr);
end;

_ParallelBinaryMatrixShapesByDimWithPrints := function(nr_rows, nr_cols,
                                                       fork_nr, nr_forks)
  local data, range, nr_data, report_gap, count, last_count, embed, G,
  extensions, out, temp, i, e;

  # If we have calculated this case before then return those results
  data := _ParallelBinaryMatrixShapesReadFile(nr_rows, nr_cols);
  if data <> fail then
    return data;
  fi;

  # In square matrix case, build up from previous square case
  if nr_rows = nr_cols then
    data := _ParallelBinaryMatrixShapesReadFile(nr_rows - 1, nr_cols - 1);
    range := Int(Floor(Float(Size(data) / nr_forks)));
    if fork_nr = nr_forks then
      range := [1 + range * (nr_forks - 1) .. Size(data)];
    else
      range := [1 + range * (fork_nr - 1) .. range * fork_nr];
    fi;
    nr_data := Size(range);
    if nr_data >= 1000 then
      report_gap := 5;
    elif nr_data >= 100 then
      report_gap := 2;
    else
      report_gap := 1;
    fi;
    count := 1;
    last_count := 1;

    embed := EmbedSetTransformation(nr_rows - 1, nr_cols - 1, nr_rows, nr_cols);
    Apply(data, d -> OnSets(d, embed));
    G := ActingBinaryMatrixGroup(nr_rows, nr_cols);
    extensions := SquareExtensions(nr_rows);
    out := [];
    for i in range do
      if count = last_count + report_gap then
        Print("At data position ", String(count), " of ", String(nr_data),
        "\n");
        last_count := count;
      fi;
      count := count + 1;
      for e in extensions do
        temp := ShallowCopy(data[i]);
        UniteSet(temp, e);
        AddSet(out, CanonicalImage(G, temp, OnSets));
      od;
    od;

  # If more rows than columns then build from a square matrix by adding rows.
  elif nr_rows > nr_cols then
    data := _ParallelBinaryMatrixShapesReadFile(nr_rows - 1, nr_cols);
    range := Int(Floor(Float(Size(data) / nr_forks)));
    if fork_nr = nr_forks then
      range := [1 + range * (nr_forks - 1) .. Size(data)];
    else
      range := [1 + range * (fork_nr - 1) .. range * fork_nr];
    fi;
    nr_data := Size(range);
    if nr_data >= 1000 then
      report_gap := 5;
    elif nr_data >= 100 then
      report_gap := 2;
    else
      report_gap := 1;
    fi;
    count := 1;
    last_count := 1;

    G := ActingBinaryMatrixGroup(nr_rows, nr_cols);
    extensions := RowExtensions(nr_rows, nr_cols);
    out := EmptyPlist(Size(extensions) * Size(data));
    for i in range do
      if count = last_count + report_gap then
        Print("At data position ", String(count), " of ", String(nr_data),
        "\n");
        last_count := count;
      fi;
      count := count + 1;
      for e in extensions do
        temp := ShallowCopy(data[i]);
        UniteSet(temp, e);
        AddSet(out, CanonicalImage(G, temp, OnSets));
      od;
    od;
  else
    return fail;
  fi;

  _ParallelBinaryMatrixShapesWriteFile(nr_rows, nr_cols, out, fork_nr);
end;

# Collating output

# prefix := "/circa/home/cr66/shapes/shapes-6x6-fork"

# Collates files call Concatenation(prefix, "f", String(n)) where n in [1 ..
# nr_forks]
_ParallelOutputCollation := function(prefix, nr_forks)
  local out, name, file, tmp, file_out, i;
  out := [];
  for i in [1 .. nr_forks] do
    name := Concatenation(prefix, String(i));
    file := IO_File(name, "r");
    if file = fail then
      return fail;
    fi;
    tmp := IO_ReadLines(file);
    Apply(tmp, a -> EvalString(Chomp(a)));
    out := Union(out, tmp);
    IO_Close(file);
  od;

  file_out := IO_File(prefix, "w");
  IO_WriteLines(file_out, out);
  IO_Close(file_out);
end;

_Parallelfoo4 := function(n, filename)
  local _AsBipartiteDigraph, _AsBooleanMat, out, maps, digraphs, groups, lookup,
  reps, last_print, G, j, pos, size, really_out, gr1, zeros, nbs, bp, p, gr2,
  sort, i, k, rep, v;

  # Convert a digraph to a bipartite digraph via the adjacency matrix with all
  # edges pointing in the same direction (so that we may act on the rows and
  # columns independently)
  _AsBipartiteDigraph := function(nbs)
    return DigraphNC(Concatenation(List(nbs, x -> x + n),
                                   List([1 .. n], x -> [])));
  end;

  # convert a digraph (of the type created by _AsBipartiteDigraph above) back
  # into a boolean mat.
  # pair[1] is a canonical bipartite digraph; pair[2] is the canonical labelling
  # that made pair[1].
  _AsBooleanMat := function(pair)
    local n, ver, out, mat, inv, i;

    n := DigraphNrVertices(pair[1]) / 2;
    ver := OnTuples([1 .. n], pair[2]);
    out := OutNeighbours(pair[1]);
    mat := EmptyPlist(n);
    inv := pair[2] ^ -1;

    for i in [1 .. n] do # consider how to calculate pair[2] ^-1 less often
      mat[i] := BlistList(ver, OnTuples(OnTuples(out[ver[i]], inv) - n,
      pair[2]));
    od;

    return MatrixNC(BooleanMatType, mat);
  end;

  out := [];
  maps := [];
  Print("Reading digraphs from file . . .\n");
  digraphs := ReadDigraphs(filename, TCodeDecoderNC);
  groups := [];
  lookup := [];
  reps   := [];
  last_print := 0;

  Print("Finding automorphism groups of digraphs . . .\n");
  for i in [1 .. Length(digraphs)] do
    if i = last_print + 1000 then
      Print("At ", i, " of ", Length(digraphs), ", found ");
      Print(Length(groups), " automorphism groups, so far\n");
      last_print := i;
    fi;
    G := AutomorphismGroup(digraphs[i]);
    j := Position(groups, G);

    if j = fail then
      Add(groups, G);
      pos       := Length(groups);
      lookup[i] := pos;
      reps[pos] := [];
      for k in [1 .. n] do
        Add(reps[pos], List(Orbits(G, Combinations([1 .. n], k), OnSets),
                            Representative));
      od;
    else
      lookup[i] := j;
    fi;
  od;
  Print(Length(groups), " automorphism groups found\n");

  size := 0;
  really_out := [];
  Print("Finding representatives . . .\n");
  last_print := 0;
  for i in [1 .. Length(digraphs)] do
    if i = last_print + 999 then
      Print("At ", i, " of ", Length(digraphs));
      Print(", found ", Length(really_out), " representatives, so far\n");
      last_print := i;
    fi;
    pos := lookup[i];
    gr1 := digraphs[i];
    zeros := Filtered([1 .. n], v -> OutDegrees(gr1)[v] = 0
                                     or InDegrees(gr1)[v] = 0);

    for j in [Maximum(Size(zeros), 1) .. n] do
      for rep in reps[pos][j] do
        if IsSubsetSet(rep, zeros) then
          # No rows or columns of zeros in adjacency matrix of gr with a loop
          # added to every vertex in rep.
          nbs := OutNeighboursMutableCopy(gr1);
          for v in rep do
            AddSet(nbs[v], v);
          od;
          bp := _AsBipartiteDigraph(nbs);
          p := BlissCanonicalLabelling(bp);
          gr2 := OnDigraphs(bp, p);

          sort := PositionSorted(out, gr2);
          if sort > size then
            Add(out, gr2);
            if BlissCanonicalDigraph(DigraphReverse(gr2)) >= gr2 then
              Add(really_out, _AsBooleanMat([gr2, p]));
            fi;
            size := size + 1;
          elif gr2 <> out[sort] then
            CopyListEntries(out, sort, 1, out, sort + 1, 1, size);
            out[sort] := gr2;
            if BlissCanonicalDigraph(DigraphReverse(gr2)) >= gr2 then
              Add(really_out, _AsBooleanMat([gr2, p]));
            fi;
            size := size + 1;
          fi;
        fi;
      od;
    od;
  od;
  return really_out;
end;

#########
# Testing
#########

shortcutss := function(nr_rows, nr_cols)
  local data, embed, G, extensions, out, pos, temp, d, e;

  # In square matrix case, build up from previous square case
  if nr_rows = nr_cols then
    data := BinaryMatrixShapesByDim(nr_rows - 1, nr_cols - 1);
    embed := EmbedSetTransformation(nr_rows - 1, nr_cols - 1, nr_rows, nr_cols);
    Apply(data, d -> OnSets(d, embed));
    G := ActingBinaryMatrixGroup(nr_rows, nr_cols);
    extensions := SquareExtensions(nr_rows);
    out := [];
    for d in data do
      pos := Position(data, d);
      for e in extensions do
        Print("\c Rows = ", nr_rows, "; Cols = ", nr_cols, "; Data pos = ", pos,
        "/", Size(data), "; Extensions pos = ", Position(extensions, e), "/",
        Size(extensions), "\r");
        temp := ShallowCopy(d);
        UniteSet(temp, e);
        AddSet(out, CanonicalImage(G, temp, OnSets));
      od;
    od;

  # If more rows than columns then build from a square matrix by adding rows.
  elif nr_rows > nr_cols then
    data := BinaryMatrixShapesByDim(nr_rows - 1, nr_cols);
    G := ActingBinaryMatrixGroup(nr_rows, nr_cols);
    extensions := RowExtensions(nr_rows, nr_cols);
    out := EmptyPlist(Size(extensions) * Size(data));
    for d in data do
      pos := Position(data, d);
      for e in extensions do
        Print("\c Rows = ", nr_rows, "; Cols = ", nr_cols, "; Data pos = ", pos,
        "/", Size(data), "; Extensions pos = ", Position(extensions, e), "/",
        Size(extensions), "\r");
        temp := ShallowCopy(d);
        UniteSet(temp, e);
        AddSet(out, CanonicalImage(G, temp, OnSets));
      od;
    od;
  fi;

  return out;
end;
