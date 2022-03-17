#!/bin/bash

echo "./deploy.sh $*" > redeploy.sh
chmod +x redeploy.sh

existing=$(docker ps -aqf name=smiles-loader)
if [ -n "$existing" ]; then
    echo "removing existing container"
    docker rm -f $existing
fi

docker run -d \
--name smiles-loader \
--restart unless-stopped \
-e ARGS="$*" \
smiles-loader
