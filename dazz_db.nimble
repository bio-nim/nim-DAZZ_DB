# Package

version       = "0.0.0"
author        = "Christopher Dunn"
description   = "Nim wrapper for DAZZ_DB (Gene Myers)"
license       = "MIT"

# Dependencies

requires "nim >= 0.17.0"

srcDir = "./src"

if not fileExists("repos/DAZZ_DB/DB.h"):
    let msg = "git submodule update --init"
    echo msg
    exec(msg)

task test, "Test dazz_db wrapper":
    withDir("tests"):
        exec("make")
