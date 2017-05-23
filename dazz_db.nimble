# Package

version       = "0.0.0"
author        = "Christopher Dunn"
description   = "Nim wrapper for DAZZ_DB (Gene Myers)"
license       = "MIT"

# Dependencies

requires "nim >= 0.17.0"

srcDir = "./src"

before install:
    echo "HELLO"

task test, "Test dazz_db wrapper":
    withDir("tests"):
        exec("nim c -r open")
