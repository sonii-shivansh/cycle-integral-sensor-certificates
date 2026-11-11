# Release and Zenodo archival guide

This repository is a preliminary **scientific software prototype**, not a peer-reviewed paper or a claim of scientific novelty.

## Before creating a release

1. Check that [LICENSE](LICENSE) is the intended license. The repository currently uses **MIT**.
2. Check the creator's public name, version, and description in [CITATION.cff](CITATION.cff). Its version field is prepared for **0.1.0**.
3. Review the proofs, limitations, and prior-art references in [CIFC_research_note.md](CIFC_research_note.md), [LIMITATIONS.md](LIMITATIONS.md), and [REFERENCES.md](REFERENCES.md).
4. Run the demonstration and regression tests:

   ~~~bash
   python -m pip install -r requirements.txt
   python -m unittest discover -s tests -v
   python cycle_integral_fault_certificate.py
   ~~~

5. **Before** creating a GitHub release, sign in at https://zenodo.org, link GitHub, use **Sync now**, and enable this exact repository in Zenodo's GitHub integration. Otherwise an existing release may not be archived automatically.

## Create the first release

On GitHub, open:
https://github.com/sonii-shivansh/cycle-integral-sensor-certificates/releases/new

- Create a new tag: **v0.1.0**
- Target branch: **main**
- Title: **CIFC v0.1.0 — Initial Research Prototype**
- Include that the research is preliminary, unreviewed, and based on synthetic results.
- Publish the release.

Zenodo should archive newly published releases from **enabled repositories**. Verify the release appears in your Zenodo GitHub integration and inspect the publication metadata and actual assigned DOI. Do **not** invent or prefill a DOI.

## After the release receives a DOI

1. Record the **version-specific DOI** in the GitHub release notes, and optionally in citation metadata (updating citation metadata after a release changes the repository, **not** the already archived release).
2. Consider adding the **concept DOI** (which resolves to the latest archived version) to the README badge and project overview.
3. In the README, accurately explain whether any DOI points to a fixed version or the family of versions.
4. Preserve historical release tags; publish corrections as a subsequent version rather than silently rewriting an archived artifact.

An archived DOI makes software citable; it does not independently validate mathematics or establish priority over earlier research.

## Sources

- Zenodo: https://help.zenodo.org/docs/github/enable-repository/
- Zenodo citation metadata: https://help.zenodo.org/docs/github/describe-software/citation-file/
- GitHub CITATION files: https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-citation-files
