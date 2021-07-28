We vendor the entirety of https://github.com/wildart/SmithNormalForm.jl because it is not, and will not (see [SmithNormalForm.jl #2](https://github.com/wildart/SmithNormalForm.jl/issues/2)) be, registered in Julia's General registry.
However, in order for Crystalline to be registered in the General Registry, we cannot depend on an unregistrered package, so we need to vendor it ourselves.

## git subtree
We can use git's subtree functionality to pull down SmithNormalForm.jl and also keep it up-to-date:

Specifically, SmithNormalForm.jl's git repo was added following the strategy in https://www.atlassian.com/git/tutorials/git-subtree, with the commands:

```sh
git remote add -f wildart-snf https://github.com/wildart/SmithNormalForm.jl.git
git subtree add --prefix .vendor/SmithNormalForm wildart-snf master --squash
```

and we can check for updates (and integrate them) via:

```sh
git fetch wildart-snf master
git subtree pull --prefix .vendor/SmithNormalForm wildart-snf master --squash
```
