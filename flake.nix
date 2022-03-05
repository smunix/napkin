{
  description = "virtual environments";

  inputs.devshell.url = "github:numtide/devshell";
  inputs.flake-utils.url = "github:numtide/flake-utils";

  outputs = { self, flake-utils, devshell, nixpkgs }:
    flake-utils.lib.eachDefaultSystem (system: {
      devShell = let
        overlay = f: _: {
          napkin = f.stdenv.mkDerivation {
            name = "napkin";
            src = self;
            buildInputs = with f; [ libGL libGLU freeglut ];
            buildPhase =
              "g++ -march=native -O3 main.cpp -lglut -lGL -pthread -o napkin";
            installPhase = ''
              mkdir -pv $out/bin
              cp napkin $out/bin/
            '';
          };
        };
        pkgs = import nixpkgs {
          inherit system;
          overlays = [ devshell.overlay overlay ];
        };
      in pkgs.devshell.mkShell {
        imports = [ (pkgs.devshell.importTOML ./devshell.toml) ];
        packages = [ pkgs.napkin ];
      };
    });
}
