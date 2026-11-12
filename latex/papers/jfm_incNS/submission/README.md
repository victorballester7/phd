# JFM submission package

Submission site: <https://mc.manuscriptcentral.com/jfm> (ScholarOne; this account is separate from your Cambridge Core account).
Rebuild everything with `./make_submission.sh`. It compiles a copy of the manuscript, so it does not clash with VimTeX's `latexmk -pvc`.

## Files to upload

| ScholarOne designation      | File                                 | JFM requirement                                    | Status |
|-----------------------------|--------------------------------------|----------------------------------------------------|--------|
| Main Document               | `main_review.pdf`                    | PDF, figures embedded, **≤ 6 MB**                  | 1.2 MB, 32 pp (raw `main.pdf` is 29 MB: do not upload that one) |
| Statement on Significance   | `statement_of_significance.txt`      | 50–100 words (JFM Papers only)                     | 92 words |
| Graphical Abstract          | `graphical_abstract.jpg`             | 1 panel, 1.2:1, JPEG ≥ 300 dpi, avoid text        | 1905×1588 px, no text (λ₂ snapshot, w=19, d=4, Re=1000) |
| Cover Letter                | `cover_letter.pdf`                   | not mandated in detail                             | has TODOs |
| Movie / Other Supp. Material| none                                 | MP4 (H.264), ≤ 50 MB each, cited in text           | optional, see below |

Do **not** upload LaTeX source or separate figure files now. They are only needed after acceptance.

## Still to do (marked `[TODO: …]` in red in the PDFs)

1. **Funding statement** (`main.tex`, back matter). Give the funders and grant numbers, and mark grants held by individual authors with their initials. If there was no specific funding, use JFM's wording: *"This research received no specific grant from any funding agency, commercial or not-for-profit sectors."* ScholarOne also asks for the funding details separately.
2. **Declaration of interests**: compulsory. Confirm the suggested wording about J.D.C.'s employment at Boeing with all co-authors. If all agree there is nothing to declare, use *"The authors report no conflict of interest."*
3. **Data availability**: strongly encouraged. Deposit the session files, meshes and post-processed data (Zenodo gives a DOI) and fill in the DOI. Alternatively, replace the sentence with *"The data that support the findings of this study are available from the corresponding author upon reasonable request."*
4. **ORCIDs**: an ORCID is **mandatory for the corresponding author** in ScholarOne. Listing them in the PDF is optional; delete that block if you prefer.
5. **Author contributions**: optional. Check the drafted sentence.
6. **Acknowledgements**: I added the Imperial RCS (CX3) citation. Add ARCHER2 or other allocations if you used them.
7. **Cover letter**: add the arXiv number if you post a preprint (allowed by Cambridge's preprint policy), and optionally suggest referees or name referees to exclude.
8. **AI-use declaration**: JFM requires you to declare generative-AI use for text or images (acknowledgements, methods or captions) and to give the tool, version, dates and how it was used. Spelling and grammar checks need no declaration. Decide whether anything in the manuscript needs this.

## ScholarOne form entries

- **Article type:** JFM Papers.
- **Associate Editor:** **Yongyun Hwang is a JFM Associate Editor and a co-author, so he cannot handle the paper.** I suggest **Denis Sipp** (ONERA: global stability and resolvent analysis, open-cavity flows). Alternatives: C. P. Caulfield (Cambridge, stability/bifurcations) or Genta Kawahara (Osaka, dynamical-systems view of transition). If you choose someone else, change `\AEname` and `\AEdear` in `cover_letter.tex`.
- **Keywords (pick ≤ 3 from the dropdown, not in the manuscript):** suggested *boundary layer stability*, *transition to turbulence*, *bifurcation*. Alternatives if those are not offered: *instability*, *separated flows*, *absolute/convective instability*.
- **Title:** Transition over surface gaps in incompressible laminar boundary layers
- **Authors / affiliations:** V. Ballester Ribó¹, J. D. Crouch², Y. Hwang¹, S. J. Sherwin¹ (¹ Dept. of Aeronautics, Imperial College London; ² The Boeing Company, Seattle)
- **Abstract:** paste from `main.tex` (about 230 words; limit 250, single paragraph, must stay on page 1 ✔).
- You must confirm that a Competing Interests statement and a Funding statement appear in the PDF.
- **Open access / funder mandates:** you will be asked whether a funder requires deposit in (Europe) PubMed Central. UKRI open-access rules apply if UKRI (EPSRC) funded the work.

## Changes made to `main.tex`

- Removed `\usepackage{showframe}`, which drew layout frames on every page.
- Added the JFM back-matter statements before the references (Acknowledgements, Funding, Declaration of interests, Data availability, ORCIDs, Author contributions), using the class's `bmhead` environment.

## Optional extras

- **Movies:** an animation of the transition sequence in figure 18 (w=19, d=4) or of the region-III breakdown (w=70, d=1.5) would fit JFM's movie format. Each movie must be referenced in the text as "movie 1", and so on.
- **Alt-text:** JFM generates alt-text automatically at proof stage unless you submit the Accessibility Descriptions form. Check it at proofs.
- **At acceptance:** you will need the LaTeX source, the individual figure files (named `Fig1.pdf`, …), a final PDF and the signed licence to publish.
