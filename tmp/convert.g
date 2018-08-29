

ReadCCRFile := function(fname)
  local file, out;
  file := IO_File(UserHomeExpand(fname), "r");
  if file = fail then
    return fail;
  fi;
  out := IO_ReadLines(file);
  IO_Close(file);
  Apply(out, a -> EvalString(Chomp(a)));
  return out;
end;

WriteJDMFile := function(fname, vals)
  local file, str, list, x;
  fname := Concatenation(fname, ".gz");
  file  := IO_CompressedFile(UserHomeExpand(fname), "w");
  if file = fail then
    return fail;
  fi;
  for list in vals do
    str := "";
    for x in list do
      Add(str, CharInt(x + 97));
    od;
    IO_WriteLine(file, str);
  od;
  IO_Close(file);
  return IO_OK;
end;

ReadJDMFile := function(fname)
  local file, str, out, pos, line, char;
  fname := Concatenation(fname, ".gz");
  file  := IO_CompressedFile(UserHomeExpand(fname), "r");
  if file = fail then
    return fail;
  fi;
  str := List(IO_ReadLines(file), Chomp);
  IO_Close(file);
  out := List([1 .. Length(str)], i -> EmptyPlist(Length(str[i])));
  pos := 0;
  for line in str do
    pos := pos + 1;
    for char in line do
      Add(out[pos], INT_CHAR(char) - 97);
    od;
  od;
  return out;
end;

ConvertCCRToJDM := function(fname)
  local vals;
  vals := ReadCCRFile(fname);
  return WriteJDMFile(fname, vals);
end;

dir := DirectoryContents(UserHomeExpand("~/git/0-simple-semigroups/shapes"));
dir := Filtered(dir, x -> StartsWith(x, "shapes"));
prefix := UserHomeExpand("~/git/0-simple-semigroups/shapes");

for fname in dir do 
  ConvertCCRToJDM(Concatenation(prefix, "/", fname));
od;
