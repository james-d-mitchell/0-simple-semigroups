###############################################################################
# A matrix with entries from a 0-group is represented as a point in 3D space,
# Rows X Cols X G u {0}.
###############################################################################
# I assume points in 3D space are represented as [x,y,z], where x runs from
# [1..dimensions[1]], y runs from [1..dimensions[2]], and z runs from
# [1..dimensions[3]]. These points may be represented in 1D space
# [1..dimensions[1] * dimensions[2] * dimensions[3]]. The representation in 1D
# space will be the best for our calculations as it allows us to use efficient
# methods involving permutation groups.
###############################################################################

# This code requires the images package
LoadPackage("images");

###############################################################################
# These two functions just map between 3D and 1D represenations of 0-group
# matrices, which is a little annoying as GAP counts from 1.
###############################################################################
###############################################################################
Flatten3DPoint := function(dimensions, point)
  return (point[1] - 1) * dimensions[2] * dimensions[3] +
    (point[2] - 1) * dimensions[3] + (point[3] - 1) + 1;
end;

Unflatten3DPoint := function(dimensions, value)
   local ret;
   ret := [];
   value := value - 1;
   ret[3] := value mod dimensions[3] + 1;
   value := value - (ret[3] - 1);
   value := value / dimensions[3];
   ret[2] := value mod dimensions[2] + 1;
   value := value - (ret[2] - 1);
   ret[1] := value / dimensions[2] + 1;
   return ret;
end;

###############################################################################
###############################################################################
# Helpers for working between binary matrix 'shapes' and 0-group matrices (both
# in set representations)
###############################################################################
###############################################################################
# Converts a 'flat' 2D point into the 3D point with the first two dimensions
# the same and 1 in the final dimension.
Unflatten2DPointIn3D := function(dimensions, value)
  local ret;
   ret := [];
   value := value - 1;
   ret[2] := value mod dimensions[2] + 1;
   value := value - (ret[2] - 1);
   value := value / dimensions[2];
   ret[1] := value mod dimensions[1] + 1;
   ret[3] := 1;
   return ret;
end;

###############################################################################
# Converts binary matrix to a group matrix over the group {0, 1}
BinaryMatrixToZeroGroupMatrix := function(shape, dimensions)
  local 3dshape, point, i;
  3dshape := [];
  for i in [1 .. dimensions[1] * dimensions[2]] do
    point := Unflatten2DPointIn3D(dimensions, i);
    if i in shape then
      point[3] := 2;
      Add(3dshape, point);
    else
      Add(3dshape, point);
    fi;
  od;
  Apply(3dshape, a -> Flatten3DPoint(dimensions, a));
  return 3dshape;
end;

###############################################################################
# Converts 3D matrix into a 'flat' 2D point
Flatten3DPointIn2D := function(dimensions, point)
  return (point[1] - 1) * dimensions[2] + point[2];
end;

###############################################################################
###############################################################################
# Go between 1D and 3D repesentations of 0-groups matrices
###############################################################################
###############################################################################
SetToZeroGroupMatrix := function(set, nr_rows, nr_cols, G)
  local 0G, mat, dim, point, x;
  0G := Concatenation([0], Elements(G));
  mat := List([1 .. nr_rows], a -> EmptyPlist(nr_cols));
  dim := [nr_rows, nr_cols, Size(G) + 1];
  for x in set do
    point := Unflatten3DPoint(dim, x);
    mat[point[1]][point[2]] := 0G[point[3]];
  od;
  return mat;
end;

ZeroGroupMatrixToSet := function(mat, nr_rows, nr_cols, G)
  local set, dim, i, j;
  set := [];
  dim := [nr_rows, nr_cols, Size(G) + 1];
  for i in [1 .. nr_rows] do
    for j in [1 .. nr_cols] do
      if mat[i][j] = 0 then
        Add(set, Flatten3DPoint(dim, [i, j, 1]));
      else
        Add(set,
          Flatten3DPoint(dim, [i, j, Position(Elements(G), mat[i][j])]));
      fi;
    od;
  od;
  return set;
end;

###############################################################################
###############################################################################
# The functions return permutations corresponding to actions on the space
# of 0-group matrices which have orbits corresponding to isomorphism classes
# of RZMS. In particular, these actions are: row and column permutations, and
# multiplication of rows and columns by non-zero group elements.
# The action which transposes matrices sends 0-group matrices to
# anti-isomorphic matrices.
###############################################################################
###############################################################################
# Given a 3D cube, dimensions dimensions, apply permutation 'perm' to dimension
# 'dim'.
#
# With dim = 1 and dim = 2 this returns a permutation which acts on 1D 0-group
# matrices by permuting the rows and columns, respectively.
# With dim = 3 and perm corresponding to an automorphism of the 0-group this
# returns a permutation on 1D reps which applies an automorphism elementwise.
ApplyPermWholeDimension := function(dimensions, dim, perm)
  local map, point, i;
    map := [];
    for i in [1 .. Product(dimensions)] do
      point := Unflatten3DPoint(dimensions, i);
      point[dim] := point[dim] ^ perm;
      map[i] := Flatten3DPoint(dimensions, point);
    od;
  return PermList(map);
end;

###############################################################################
# Given a 3D cube, dimensions dimensions, apply permutation 'perm' to dimension
# 'dim' But only to elements which have value 'fixval' in dimension 'fixdim'.
#
# Returns a permutation of 1D 0-group matrices corresponding to multiplying a
# row (fixdim = 1) or column (fixdim = 2) when:
#
# (i)  dim = 3, and
# (ii) perm is a permutation of [1..dimensions[3]] corresponding to the action
# (by left or right multiplication) of a non-zero element of G_0.
#
# (That is to say, the elements of G_0 correspond to [1 .. Size(G_0)] and perm
# sends a -> a * b for all a in G_0 and some b in G_0.)
ApplyPermSingleAssignDimension  := function(dimensions, dim, perm, fixdim,
                                            fixval)
  local map, point, i;
  map := [];
  for i in [1 .. Product(dimensions)] do
    point := Unflatten3DPoint(dimensions, i);
    if point[fixdim] = fixval then
      point[dim] := point[dim] ^ perm;
    fi;
    map[i] := Flatten3DPoint(dimensions, point);
  od;
  return PermList(map);
end;

###############################################################################
# Returns a permutation of 1D 0-group matrices with dimensions dimensions
# corresponding to transposition of matrices.
TranspositionPerm := function(dimensions)
  local map, point, i;
  map := [];
  for i in [1 .. Product(dimensions)] do
    point := Unflatten3DPoint(dimensions, i);
    map[i] := Flatten3DPoint(dimensions, [point[2], point[1], point[3]]);
  od;
  return PermList(map);
end;

###############################################################################
# This function uses ApplyPermWholeDimension and ApplyPermSingleAssignDimension
# to construct the group which has orbits corresponding to distinct isomorphism
# classes of RZMS with a certain 'shape'. Having the same 'shape' means that
# RZMS are created from matrices which have zeros in the same locations.
RZMSMatrixIsomorphismGroup := function(shape, nr_rows, nr_cols, G)
  local dim, S, rows, cols, 3dshape, H, gens, elms, rmlt, grswaps, lmlt,
    gcswaps, auto, g;

  dim := [nr_rows, nr_cols, Size(G) + 1];
  # Row Swaps
  S := SymmetricGroup(nr_rows);
  rows := List(GeneratorsOfGroup(S), x -> ApplyPermWholeDimension(dim, 1, x));

  # Col swaps
  S := SymmetricGroup(nr_cols);
  cols := List(GeneratorsOfGroup(S), x -> ApplyPermWholeDimension(dim, 2, x));

  # We only want swaps which fix the locations of the zeros.
  3dshape := BinaryMatrixToZeroGroupMatrix(shape, dim);
  H := Stabilizer(Group(Flat([cols, rows])), 3dshape, OnSets);

  gens := GeneratorsOfGroup(G);
  elms := ShallowCopy(Elements(G));

  # Apply g to each row (right multiplication):
  rmlt := List(gens, g -> PermList(Concatenation([1],
          1 + List(elms, e -> Position(elms, e * g)))));
  grswaps := List([1 .. dim[1]], r -> List(rmlt, g ->
  ApplyPermSingleAssignDimension(dim, 3, g, 1, r)));

  # Apply g to each col (left multiplication by inverse):
  lmlt := List(gens, g -> PermList(Concatenation([1],
          1 + List(elms, e -> Position(elms, g ^ -1 * e)))));
  gcswaps := List([1 .. dim[2]], r -> List(lmlt, g ->
  ApplyPermSingleAssignDimension(dim, 3, g, 2, r)));

  # Automorphisms of G
  S := AutomorphismGroup(G);
  auto := List(GeneratorsOfGroup(S), x -> List(Elements(G), a ->
          Position(Elements(G), a ^ x)));
  Apply(auto, a -> PermList(Concatenation([1], a + 1)));
  auto := List(auto, x -> ApplyPermWholeDimension(dim, 3, x));

  # The RZMS matrix isomorphism group
  g := Group(Flat([GeneratorsOfGroup(H), grswaps, gcswaps, auto]));

  return g;
end;

###############################################################################
###############################################################################
# This function finds some of the possible dimensions of matrices which define
# RZMS of order k. I will find [x, y, z] such that x * y * z + 1 = k,
# however if [x, y, z] and [y, x, z] are distinct then it will only return one
# of these because the isomorphism classes of RZMS with dimensions [x, y, z]
# will be in 1-1 correspondence (by anti-isomorphism) with the classes of RZMS
# which have dimensions [y, x, z]. It is easy to go between these classes as we
# only need to transpose the matrices so there is no point calculating both
# cases by finding orbit representatives of a potentially very large
# permutation group.
#
# If anti_iso = true then we look up to anti isomrophism, otherwise just up to
# isomorphism.
###############################################################################
###############################################################################
FindRZMSTripleParametersByOrder := function(k, anti_iso)
  local out, d, e;
  out := [];
  for d in DivisorsInt(k - 1) do
    # |I| * |J| = k/d
    for e in DivisorsInt((k - 1) / d) do
      # If anti = true then We only want one of (|I|,|J|) = (a,b), (b,a)
      if not anti_iso or e <= (k - 1) / d / e then
        Add(out, [(k - 1) / d / e, e, d]);
      fi;
    od;
  od;
  return out;
end;

###############################################################################
###############################################################################
# Methods for enumerating RZMS matrices
###############################################################################
###############################################################################
# Given dimensions and a group and a shape, returns a representative for each
# isomorphism class of RZMS defined by a matrix with these dimensions over the
# 0-group G u {0} with 0's precisely the locations specified by the shape.
# When nr_rows is not equal to nr_cols a representative of each
# anti-isomorphism class is returned instead.
#
# We reduce the number of matrices to apply CanonicalImage to by only applying
# to certain normal matrices, which will still cover every isomorphism class.
RZMSMatricesByShape := function(nr_rows, nr_cols, G, shape)
  local IG, dim, space, is_norm_row, is_norm_col, point, out, iter, mat, i;

  # This is slow when one of nr_rows or nr_cols is large (even if the other is
  # small) so we have RZMSMatricesByShapeEasyCase for some of these situations.
  IG := RZMSMatrixIsomorphismGroup(shape, nr_rows, nr_cols, G);
  dim := [nr_rows, nr_cols, Size(G) + 1];
  space := [];
  is_norm_row := List([1 .. nr_rows], a -> false);
  is_norm_col := List([1 .. nr_cols], a -> false);
  for i in [1 .. dim[1] * dim[2]] do
    point := Unflatten2DPointIn3D(dim, i);
    if i in shape then
      if not is_norm_row[point[1]] then
        Add(space, [Flatten3DPoint(dim, [point[1], point[2], 2])]);
        is_norm_row[point[1]] := true;
      elif not is_norm_col[point[2]] then
        Add(space, [Flatten3DPoint(dim, [point[1], point[2], 2])]);
        is_norm_col[point[2]] := true;
      else
        Add(space, List([1 .. Size(G)], g ->
          Flatten3DPoint(dim, [point[1], point[2], g + 1])));
      fi;
    else
      Add(space, [Flatten3DPoint(dim, [point[1], point[2], 1])]);
    fi;
  od;

  out := Set([]);
  iter := IteratorOfCartesianProduct(space);
  while not IsDoneIterator(iter) do
    mat := NextIterator(iter);
    AddSet(out, CanonicalImage(IG, mat, OnSets));
  od;

  return out;
end;

###############################################################################
# Given a trivial case (nr_rows = 1 or nr_cols = 1 or Size(G) = 1) this returns
# the (anti-)isomorphism classes of RZMS with the given dimensions. The method
# RZMSMatricesByShape would be very ineffective for these cases.
RZMSMatricesByShapeEasyCase := function(nr_rows, nr_cols, G, shape)
  local out, dim, pos, point, i;

  out := [];
  dim := [nr_rows, nr_cols, Size(G) + 1];
  if dim[1] = 1 or dim[2] = 1 then
    pos := Position(Elements(G), One(G)) + 1;
    return [List([1 .. nr_rows], i -> pos + (i - 1) * (Size(G) + 1))];
  elif dim[3] = 2 then
    for i in [1 .. dim[1] * dim[2]] do
      point := Unflatten2DPointIn3D(dim, i);
      if i in shape then
        Add(out, Flatten3DPoint(dim, [point[1], point[2], 2]));
      else
        Add(out, Flatten3DPoint(dim, [point[1], point[2], 1]));
      fi;
    od;
    return [out];
  fi;
  ErrorNoReturn("RZMSMatricesByShapeEasyCase: only for cases with 1 row, 1 ",
                "column or trivial group");
end;

###############################################################################
# Given dimensions and a group, returns a representative for each isomorphism
# class of RZMS defined by a matrix with these dimensions over the 0-group G u
# {0}. When nr_rows is not equal to nr_cols a representative of each
# anti-isomorphism class is returned instead
RZMSMatricesByParameters := function(nr_rows, nr_cols, G)
  local pos, out, i, shapes, shape;

  if nr_cols = 1 then
    pos := Position(Elements(G), One(G)) + 1;
    return [List([1 .. nr_rows], i -> pos + (i - 1) * (Size(G) + 1))];
  fi;

  # The m x n case is deduced from the n x m case.
  if nr_rows < nr_cols then
    out := RZMSMatricesByParameters(nr_cols, nr_rows, G);
    Apply(out, a -> Unflatten3DPoint([nr_rows, nr_cols, Size(G) + 1], a));
    for i in out do
      i := [i[2], i[1], i[3]];
    od;
    Apply(out, a -> Flatten3DPoint([nr_rows, nr_cols, Size(G) + 1], a));
    return out;
  fi;

  # Get shapes - if 3rd parameter is true then using precalculated binary mats.
  shapes := BinaryMatrixShapesByDim(nr_rows, nr_cols, true);

  # Find by shape
  out := [];
  for shape in shapes do
    if Size(G) = 1 then
      Append(out, RZMSMatricesByShapeEasyCase(nr_rows, nr_cols, G, shape));
    else
      Append(out, RZMSMatricesByShape(nr_rows, nr_cols, G, shape));
    fi;
  od;

  return out;
end;

###############################################################################
# Given an order, returns a representative matrix for each isomorphism class of
# RZMS, with that order, when nr_rows equals nr_cols and a representative of
# each anti-isomorphism class otherwise.
RZMSMatrices := function(order)
  local out, parameters, H, mats, p, G;
  out := [];
  parameters := FindRZMSTripleParametersByOrder(order, true);
  for p in parameters do
    for G in AllSmallGroups(p[3]) do
      H := Image(IsomorphismPermGroup(G));
      mats := RZMSMatricesByParameters(p[1], p[2], H);
      Apply(mats, mat -> SetToZeroGroupMatrix(mat, p[1], p[2], H));
      Add(out, Concatenation([H], mats));
    od;
  od;
  return out;
end;

# Enumerates and returns representatives of each isomorphism class as a RZMS
RZMS := function(order)
  local out, parameters, H, mats, p, G;
  out := [];
  parameters := FindRZMSTripleParametersByOrder(order, true);
  for p in parameters do
    for G in AllSmallGroups(p[3]) do
      H := Image(IsomorphismPermGroup(G));
      mats := RZMSMatricesByParameters(p[1], p[2], H);
      Apply(mats, mat -> SetToZeroGroupMatrix(mat, p[1], p[2], H));
      Apply(mats, mat -> ReesZeroMatrixSemigroup(H, mat));
      Add(out, mats);
    od;
  od;
  return Concatenation(out);
end;

# Just for interest (have better method elsewhere) returns just the rees matrix
# semigroups. They correspond to the RZMS which have no zeros in the defining
# matrix and are obtained by removing the zero element from these.
RMSMatrices := function(order)
  local out, parameters, H, shape, mats, p, G;
  out := [];
  parameters := FindRZMSTripleParametersByOrder(order, true);
  for p in parameters do
    for G in AllSmallGroups(p[3]) do
      H := Image(IsomorphismPermGroup(G));
      shape := [1 .. p[1] * p[2]];
      if p[1] = 1 or p[2] = 1 or Size(H) = 1 then
        mats := RZMSMatricesByShapeEasyCase(p[1], p[2], H, shape);
      else
        mats := RZMSMatricesByShape(p[1], p[2], H, shape);
      fi;
      Apply(mats, mat -> SetToZeroGroupMatrix(mat, p[1], p[2], H));
      Add(out, [H, mats]);
    od;
  od;
  return out;
end;

RZMSMain := function(arg)
  local order, row, col, group_order, group, anti_iso, param, out, grps, mats, p, G;

  if Length(arg) > 5 or Length(arg) = 0 then
    ErrorNoReturn("RZMS: there should be 1, 2, 3, 4 or 5 arguments");
  fi;

  if Length(arg) >= 1 then
    if IsInt(arg[1]) and arg[1] > 1 then
      order := arg[1];
    else
      ErrorNoReturn("RZMS: first argument should be a positive integer or 0,");
    fi;
  fi;

  if Length(arg) >= 2 then
    if IsInt(arg[2]) and arg[2] >= 0 then
      row := arg[2];
    else
      ErrorNoReturn("RZMS: second argument should be a positive integer ",
                    "or 0,");
    fi;
  else
    row := 0;
  fi;

  if Length(arg) >= 3 then
    if IsInt(arg[3]) and arg[3] >= 0 then
      col := arg[3];
    else
      ErrorNoReturn("RZMS: third argument should be a positive integer or 0,");
    fi;
  else
    col := 0;
  fi;

  if Length(arg) >= 4 then
    if IsInt(arg[4]) and arg[4] >= 0 then
      group := arg[4];
    elif IsPermGroup(arg[4]) then
      group := arg[4];
    else
      ErrorNoReturn("RZMS: fourth argument should be a positive integer, or ",
                    "0, or a perm group,");
    fi;
  else
    group := 0;
  fi;
  
  if Length(arg) = 5 then
    if IsBool(arg[5]) then 
      anti_iso := arg[5]; # If true then we enumerate up to anti isomorphism
    else
      ErrorNoReturn("RZMS: fifth argument should true or false,");
    fi;
  else
    anti_iso := true; # Default setting
  fi;

  if IsInt(group) then
    param := SEMIGROUPS.RZMSParameters(order, row, col, group, anti_iso);
  else
    param := SEMIGROUPS.RZMSParameters(order, row, col, Size(group), anti_iso);
  fi;

  out := [];
  for p in param do
    if IsInt(group) then
      grps := List(AllSmallGroups(p[3]), G -> Image(IsomorphismPermGroup(G)));;
    else
      grps := [group];
    fi;
    for G in grps do
      mats := RZMSMatricesByParameters(p[1], p[2], G);
      Apply(mats, mat -> SetToZeroGroupMatrix(mat, p[1], p[2], G));
      Apply(mats, mat -> ReesZeroMatrixSemigroup(G, mat));
      Add(out, mats);
    od;
  od;
  return Concatenation(out);
end;

SEMIGROUPS.RZMSParameters := function(order, om, on, oG, anti_iso)
  local out, nr_rows, ord_grp, nr_cols;
  out := [];
  for ord_grp in DivisorsInt(order - 1) do
    if oG = 0 or oG = ord_grp then
      for nr_cols in DivisorsInt((order - 1) / ord_grp) do
        if on = 0 or on = nr_cols then
          nr_rows := ((order - 1) / ord_grp) / nr_cols;
          if (not anti_iso or nr_cols <= nr_rows)
              and (om = 0 or om = nr_rows) then
            Add(out, [nr_rows, nr_cols, ord_grp]);
          fi;
        fi;
      od;
    fi;
  od;
  return out;
end;
