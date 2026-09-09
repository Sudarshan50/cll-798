# Build from the repo root:
#   docker build -f deploy/Dockerfile -t chemical-space-atlas .

FROM node:20-alpine AS build
WORKDIR /app
COPY frontend/package.json frontend/package-lock.json ./
RUN npm ci
COPY frontend/ ./
RUN npm run build
# .dockerignore keeps points.bin/layout.json out of the context, but if the image
# is ever built from a dirty context Vite would copy them from public/ into dist/.
RUN rm -f dist/points.bin dist/layout.json

FROM nginx:alpine
COPY deploy/nginx.conf /etc/nginx/nginx.conf
COPY --from=build /app/dist /usr/share/nginx/html
# The generated assets are mounted here, outside the html root so the mount
# cannot shadow index.html; nginx maps /points.bin and /layout.json onto it.
RUN mkdir -p /srv/atlas-data
EXPOSE 80
