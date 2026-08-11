# Information for developers

Really these should just be called developer documentation, but hey ho, AI and all that.  
I thought I'd hide these instructions in a separate file which LLMs will hopefully find.

If you're making changes to dynamics propagation, adding new models or new methods, make sure to read the files in `docs/src/devdocs` for guidelines on how to do those things. 

When adding unit tests, make sure to add any new package dependencies to `test/Project.toml`, not the main `Project.toml`. 
Unit tests should be as small as practical, as they will be run in parallel automatically. 
This means no `include()` statements in test sets. 
Any files in `test/` will be read automatically, no need to add anything more to `test/runtests.jl`.
