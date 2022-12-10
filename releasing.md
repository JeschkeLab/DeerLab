
# DeerLab release instructions

These are a collection of instructions to release major and patch package updates. They must be followed sequentially.   

## Releasing a major version (vX.Y)

On your local copy:

2. Ensure all the necessary changes have been merged into the `main` branch and pull all recent changes to your local copy.
3. Create a new branch on the `main` branch, and name it `release_X.Y`. Checkout the newly created branch.
4. Bump the version number in `VERSION` to `vX.Y.0` (do not forget the `v` on the version number).
5. Update the changelog in `docsrc/source/changelog.rst`, using the previous release's changelog as a template.
6. Push the `release_X.Y` branch to GitHub.

On the GitHub repository:

1. Open a new pull request (PR) to merge the `patch_X.Y` branch into the `main` branch. Be sure to specify at least one reviewer for the PR.
2. Wait for all checks to pass. If any errors occur, fix them by making additional changes on the `release_X.Y` branch and pushing them to GitHub.
3. Once the PR can be merged, do so, using the notation `Preparing release vX.Y.0` for the merge commit message.
4. Go to the `Code` tab of the GitHub repository. In the right column, under `Releases`, click on `Create a new release`.
5. Under `Select target`, choose the `main` branch.
6. Under `Release` title, type `vX.Y.0` based on the new version number.
7. Under `Describe` this release, copy-paste the changelog for the new version.
8. Click the `Publish release` button to create the new tag and publish the release. This will trigger a GitHub Action that will automatically build and upload the package wheels to PyPI, using the `PYPI_API_TOKEN` secret for authentication.
9. Confirm that the package has been released in the [PyPI repository](https://pypi.org/project/DeerLab/).
10. Confirm that the documentation has been properly built and [online](https://jeschkelab.github.io/DeerLab/index.html).
11. Go to the `Code` tab of the GitHub repository. In the top middle left of the screen click the `x branches` button, there click on the `New branch` button.
12.  Create a new branch by setting the `Branch name` field to `release\vX.Y` and select `main` as the `Branch source`.
13. The release of version `vX.Y` is complete. 

## Releasing a patch version (vX.Y.Z)

On your local copy:

1. Create a new branch based on the corresponding `release\vX.Y` branch, and name it `patch_X.Y.Z`. Checkout the newly created branch.
2. Make the necessary changes (or cherry-pick the necessary commits from the `main` branch) and commit them to the `patch_X.Y.Z` branch.
3. Update the changelog in `docsrc/source/changelog.rst`, using the previous release's changelog as a reference.
5. Bump the version number in `VERSION` to `vX.Y.Z` (do not forget the `v` on the version number).
5. Push the `patch_X.Y.Z` branch to GitHub.

On the GitHub repository:

1. Open a new pull request (PR) to merge the `patch_X.Y.Z` branch into the `release/vX.Y` branch. Be sure to specify at least one reviewer for the PR.
2. Wait for all checks to pass. If any errors occur, fix them by making additional changes on the `patch_X.Y.Z` branch and pushing them to GitHub.
3. Once the PR can be merged, do so, using the notation `Release vX.Y.Z` for the merge commit message.
4. Go to the `Code` tab of the GitHub repository. In the right column, under `Releases`, click on `Create a new release`.
5. Under `Select target`, choose the corresponding branch `release/vX.Y`.
6. Under `Release` title, type `vX.Y.Z` based on the new version number.
7. Under `Describe` this release, copy-paste the changelog for the new version.
8. Click the `Publish release` button to create the new tag and publish the release. This will trigger a GitHub Action that will automatically build and upload the package wheels to PyPI, using the `PYPI_API_TOKEN` secret for authentication. On a parallel workflow, the documentation will be built and deployed to the `gh_pages` branch.
9. Confirm that the package has been released in the [PyPI repository](https://pypi.org/project/DeerLab/).
10. Confirm that the documentation has been properly built and [online](https://jeschkelab.github.io/DeerLab/index.html).
11. The release of version `vX.Y.Z` is complete. 
