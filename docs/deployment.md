# Deployment status

This build is designed for local single-user exploration. The API binds to `127.0.0.1` by default and serves static assets from the repository. It has no production authentication, PostgreSQL, Redis, container deployment, background job isolation, request throttling, or production-grade artifact/object storage.

Keep the API bound to loopback unless network exposure, authentication, and access controls are implemented and reviewed. Genomic files and analysis outputs can be large; the current filesystem store has no retention policy. No Docker Compose configuration or hosted deployment target is provided or verified.
