#!/usr/bin/env bash
# Run as jay in a Linux terminal; sudo will ask for the password.
set -euo pipefail
device=/dev/disk/by-uuid/5AB8A669B8A642FD
target=/media/jay/新加卷
if ! mountpoint -q "$target"; then
    if findmnt -rn -S "$device" >/dev/null; then
        echo 'The volume is already mounted elsewhere; inspect findmnt before proceeding.' >&2
        findmnt -S "$device"
        exit 1
    fi
    sudo mkdir -p "$target"
    sudo mount -t ntfs-3g -o uid=1000,gid=1000 "$device" "$target"
fi
findmnt --mountpoint "$target"
ls -ld "$target/nektar-v2026" "$target/nektar-v2026-runs"
