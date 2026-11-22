# Contributing

Bug reports, documentation improvements, and code contributions are welcome. For larger changes, open an issue first to agree on scope.

## Development setup

Requirements: Python 3.10 or newer and Node.js 20 or newer. The demo and test suite do not require Python packages or network access.

```sh
git clone https://github.com/ranajahanzaib/gene-content-atlas.git
cd gene-content-atlas
make setup
make test
make build
```

To run the app locally, use `make dev` and open the URL printed by the server. Press Ctrl+C to stop it. Real-data acquisition and Panaroo analysis require the external tools described in the README.

## Changes and pull requests

- Keep changes focused and describe user-visible behavior in the pull request.
- Add or update tests for behavior changes.
- Run `make test`, `make build`, and `node --check apps/web/app.js` before submitting.
- Document scientific assumptions, data sources, and limitations when analysis behavior changes.
- Use clear [Conventional Commit](https://www.conventionalcommits.org/) messages.
- Do not commit credentials, `.env` files, downloaded genome packages, or generated analysis directories.

Pull requests should explain the problem, the approach, and any validation that was not possible locally.
