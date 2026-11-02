# Publishing and figure notes

## Documentation stack, from sources read on 2026-11-02
No Reddit or Stack Exchange threads could be read, so this rests on project pages, issue trackers and comparison blogs. Several of those blogs sell competing products.

- Material for MkDocs entered maintenance mode in November 2025, with a pledge of about 12 months of fixes. Its successor, Zensical, reads mkdocs.yml, but reports of plugin parity are mixed.
- Docusaurus and Starlight suit JavaScript teams and Markdown or MDX sources. Neither solves LaTeX input.
- Wiki.js and BookStack are databases for team wikis. Neither suits one LaTeX source under version control.
- Pandoc is fast but lossy, tex4ht and LaTeXML are more faithful but need package support. Our test matches: only Pandoc finished both real documents.

## Figures
- research_helpers.carto: style, equal-area maps, scale bar, Crameri colour maps, Okabe-Ito categories.
- research_helpers.autotex: numbers in LaTeX come from the analysis.
- Look at every figure after building it. Three errors were found only by looking.
