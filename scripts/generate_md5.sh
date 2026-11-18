#!/usr/bin/env bash
find $1 -type f ! -name 'MANIFEST.md5' -print0 \
  | sort -z \
  | xargs -0 md5sum > $1"/MANIFEST.md5"