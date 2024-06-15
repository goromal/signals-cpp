let
  pkgs = import (fetchTarball
    ("https://github.com/goromal/anixpkgs/archive/refs/tags/v5.15.0.tar.gz"))
    { };
in with pkgs;
clangStdenv.mkDerivation {
  name = "PACKAGENAME";
  version = "0.0.0";
  src = lib.cleanSource ./.;
  nativeBuildInputs = [ cmake ];
  buildInputs = [
    # TODO ADD deps
  ];
  preConfigure = ''
    cmakeFlags="$cmakeFlags --no-warn-unused-cli"
  '';
}
