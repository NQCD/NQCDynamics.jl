# Implementing a new model

The simplest way to understand the `NonadiabaticModels` interface is to walk through
an example implementation.

To add a new model the first step is to select the abstract type that it falls under.
After this, the new type should be created (see above for example)
and the following functions should be implemented:

As is standard in Julia, only the types shown in the leaves of the tree are concrete
types, the branches are abstract types that group similar models together.
These abstract types denote the quantities provided by each model and the functions
that are implemented.
