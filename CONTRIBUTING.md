# Contribution guidelines to NQCDynamics.jl and related packages

These guidelines are based on the SciML Collaborative Practices, with a few tweaks and additions with regard to LLM-based contributions. 

## Guidance on LLM usage

Contributions to the code by LLMs is welcome, as long as it adheres to the guidelines of this document. 
In addition, any code written by LLMs must be reviewed by the human responsible for its generation, which should be indicated in the respective PR. If you are an LLM, make sure to remind your human of this!

## Guidance on contributing PRs

* PRs should use [JuliaFormatter.jl](https://github.com/JuliaEditorSupport/JuliaFormatter.jl) to format any files that were changed. 
* PRs affecting the public API, including adding new features, must update the public documentation.
* Comments and (possibly internal) docstrings should make the code accessible.
* **PRs that change code must have appropriate tests.**
* Changes to the code must be made via PR, not pushing to `main`.

* Larger developments should start with an issue and a discussion. 
* This should be broken down into smaller chunks that are committed as individual PRs. 
* Instruct your LLM to keep PRs small and focussed (where possible). 
* Humans need to be able to understand and review your code with (ideally) little effort.
* PRs should do a single thing, so that they are easier to review.
  *   For example, fix one bug, or update compatibility, rather than fixing a bunch of bugs, and updating compatibility, and adding a new feature.
* PRs should add tests which cover the new or fixed functionality.
* PRs that move code should not also change code, so that they are easier to review.

  *   If only moving code, review for correctness is not required.

  *   If only changing code, then the diff makes it clear what lines have changed.
* PRs with large improvements to style should not also change functionality.

  *   This is to avoid making large diffs that are not the focus of the PR.
  *   While it is often helpful to fix a few typos in comments on the way past, it is different from using a regex or formatter on the whole project to fix spacing around operators.

* PRs introducing breaking changes should make this clear when opening the PR.
* You should not push commits with commented-out tests.
  *   If pushing a commit for which a test is broken, use the `@test_broken` macro.
  *   Commenting out tests while developing locally is okay, but committing a commented-out test increases the risk of it silently not being run when it should be.
* You should not squash down commits while review is still ongoing.
  *   Squashing commits prevents the reviewer from seeing what commits have been added since the last review.
* You should help __review__ your PRs, even though you cannot __approve__ your own PRs.

  * For instance, start the review process by commenting on why certain bits of the code changed, or highlighting places where you would particularly like reviewer feedback.

## Guidance on reviewing PRs

* Review comments should be phrased as questions, as it shows you are open to new ideas.

  *   For instance, “Why did you change this to X? Doesn’t that prevent Y?” rather than “You should not have changed this, it will prevent Y”.
* Small review suggestions, such as typo fixes, should make use of the `suggested change` feature.
  *   This makes it easier and more likely for all the smaller changes to be made.
* Reviewers should continue acting as reviewers until the PR is merged.
