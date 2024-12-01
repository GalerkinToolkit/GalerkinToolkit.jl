window.BENCHMARK_DATA = {
  "lastUpdate": 1733094686762,
  "repoUrl": "https://github.com/GalerkinToolkit/GalerkinToolkit.jl",
  "entries": {
    "Julia benchmark result": [
      {
        "commit": {
          "author": {
            "name": "GalerkinToolkit",
            "username": "GalerkinToolkit"
          },
          "committer": {
            "name": "GalerkinToolkit",
            "username": "GalerkinToolkit"
          },
          "id": "2c11e638a52f19e5541907eef025e7419d51e430",
          "message": "Add benchmarking suite for the handwritten poisson example",
          "timestamp": "2024-11-29T17:27:25Z",
          "url": "https://github.com/GalerkinToolkit/GalerkinToolkit.jl/pull/147/commits/2c11e638a52f19e5541907eef025e7419d51e430"
        },
        "date": 1733076786561,
        "tool": "julia",
        "benches": [
          {
            "name": "poisson/2",
            "value": 68530103,
            "unit": "ns",
            "extra": "gctime=0\nmemory=25300648\nallocs=472685\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "poisson/3",
            "value": 71401509,
            "unit": "ns",
            "extra": "gctime=0\nmemory=27994064\nallocs=532188\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "name": "GalerkinToolkit",
            "username": "GalerkinToolkit"
          },
          "committer": {
            "name": "GalerkinToolkit",
            "username": "GalerkinToolkit"
          },
          "id": "7fa1c6af93ef1e62702ab1eec107346d86eb13c4",
          "message": "Add benchmarking suite for the handwritten poisson example",
          "timestamp": "2024-11-29T17:27:25Z",
          "url": "https://github.com/GalerkinToolkit/GalerkinToolkit.jl/pull/147/commits/7fa1c6af93ef1e62702ab1eec107346d86eb13c4"
        },
        "date": 1733082337633,
        "tool": "julia",
        "benches": [
          {
            "name": "poisson/n=150",
            "value": 15012113152,
            "unit": "ns",
            "extra": "gctime=1241967732\nmemory=12825369344\nallocs=309749069\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "name": "GalerkinToolkit",
            "username": "GalerkinToolkit"
          },
          "committer": {
            "name": "GalerkinToolkit",
            "username": "GalerkinToolkit"
          },
          "id": "3f1546f904df417deef4673fbb738820bf08eb0f",
          "message": "Add benchmarking suite for the handwritten poisson example",
          "timestamp": "2024-11-29T17:27:25Z",
          "url": "https://github.com/GalerkinToolkit/GalerkinToolkit.jl/pull/147/commits/3f1546f904df417deef4673fbb738820bf08eb0f"
        },
        "date": 1733090461889,
        "tool": "julia",
        "benches": [
          {
            "name": "poisson-hand/n=10",
            "value": 7345564655,
            "unit": "ns",
            "extra": "gctime=448727698\nmemory=7481652776\nallocs=169412438\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "name": "GalerkinToolkit",
            "username": "GalerkinToolkit"
          },
          "committer": {
            "name": "GalerkinToolkit",
            "username": "GalerkinToolkit"
          },
          "id": "7000bc960bbe7d9b843985241169eccffe9d1c57",
          "message": "Add benchmarking suite with trackable performance",
          "timestamp": "2024-11-29T17:27:25Z",
          "url": "https://github.com/GalerkinToolkit/GalerkinToolkit.jl/pull/147/commits/7000bc960bbe7d9b843985241169eccffe9d1c57"
        },
        "date": 1733092863998,
        "tool": "julia",
        "benches": [
          {
            "name": "poisson-hand/n=10",
            "value": 7188105780,
            "unit": "ns",
            "extra": "gctime=388828865\nmemory=7481652776\nallocs=169412438\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "name": "GalerkinToolkit",
            "username": "GalerkinToolkit"
          },
          "committer": {
            "name": "GalerkinToolkit",
            "username": "GalerkinToolkit"
          },
          "id": "5942af9ef3a848c4d5d2862e0299cee180315f8b",
          "message": "Add benchmarking suite with trackable performance",
          "timestamp": "2024-11-29T17:27:25Z",
          "url": "https://github.com/GalerkinToolkit/GalerkinToolkit.jl/pull/147/commits/5942af9ef3a848c4d5d2862e0299cee180315f8b"
        },
        "date": 1733094686106,
        "tool": "julia",
        "benches": [
          {
            "name": "poisson-hand/n=10",
            "value": 7447960011,
            "unit": "ns",
            "extra": "gctime=430827246\nmemory=7481652776\nallocs=169412438\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      }
    ]
  }
}