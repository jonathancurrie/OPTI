# OPTI Toolbox documentation

This is the Docusaurus 3 documentation site for
[OPTI Toolbox](https://github.com/jonathancurrie/OPTI). It is configured for
GitHub Pages at `https://jonathancurrie.github.io/OPTI/`.

## Local development

From the repository root:

```powershell
Set-Location Documentation
pnpm install --frozen-lockfile
pnpm start
```

Create the production build with:

```powershell
pnpm build
```

## Migration boundary

The documentation was generated from the current public text of the legacy
PmWiki and its dedicated public image directory. PmWiki configuration,
runtime code, deleted revisions, comments, credentials, logs, and private
metadata are intentionally absent.

Legacy `/pmwiki.php/Group/Page` and `/index.php/Group/Page` routes are
generated as permanent client redirects from `migration/manifest.json`.
