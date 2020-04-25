#!/bin/sh

MIRROR_REPO=${MIRROR_REPO:-"origin"}
UPSTREAM_REMOTE=${UPSTREAM_REMOTE:-"upstream_remote"}
UPSTREAM_REMOTE_NAME=${UPSTREAM_REMOTE_NAME:-"upstream"}
UPSTREAM_REMOTE_URL=${UPSTREAM_REMOTE_URL:-"https://gitlab.erc.monash.edu.au/sembucket/semtex.git"}

echo "Pulling $UPSTREAM_REMOTE_URL"
git remote add "$UPSTREAM_REMOTE" "$UPSTREAM_REMOTE_URL"
git fetch "$UPSTREAM_REMOTE"

for branch in $(git ls-remote $UPSTREAM_REMOTE | grep -oP '(?<=refs/heads/).*'); do
    echo "Fetching branch '$branch' and setting upstream to '$MIRROR_REPO/$UPSTREAM_REMOTE_NAME/$branch'"
    git fetch "$UPSTREAM_REMOTE" "$branch:$UPSTREAM_REMOTE_NAME/$branch"
    git push --set-upstream "$MIRROR_REPO" "$UPSTREAM_REMOTE_NAME/$branch:$UPSTREAM_REMOTE_NAME/$branch"
done
