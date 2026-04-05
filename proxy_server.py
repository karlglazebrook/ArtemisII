#!/usr/bin/env python3
"""
Minimal HTTP server for Artemis II tracker.
Serves local files AND proxies JPL Horizons API requests to bypass browser CORS.
Usage: python3 proxy_server.py   (runs on port 8765)
"""
import http.server, urllib.request, urllib.parse, os, sys

PORT = 8765
HORIZONS_BASE = 'https://ssd.jpl.nasa.gov/api/horizons.api'

class Handler(http.server.SimpleHTTPRequestHandler):
    def do_GET(self):
        if self.path.startswith('/horizons-proxy?'):
            qs = self.path[len('/horizons-proxy?'):]
            url = f'{HORIZONS_BASE}?{qs}'
            try:
                req = urllib.request.Request(url, headers={'User-Agent': 'ArtemisTracker/1.0'})
                with urllib.request.urlopen(req, timeout=30) as r:
                    data = r.read()
                self.send_response(200)
                self.send_header('Content-Type', 'text/plain; charset=utf-8')
                self.send_header('Access-Control-Allow-Origin', '*')
                self.send_header('Cache-Control', 'no-store')
                self.end_headers()
                self.wfile.write(data)
            except Exception as e:
                self.send_response(502)
                self.send_header('Content-Type', 'text/plain')
                self.send_header('Access-Control-Allow-Origin', '*')
                self.end_headers()
                self.wfile.write(f'Proxy error: {e}'.encode())
        else:
            super().do_GET()

    def log_message(self, fmt, *args):
        # Suppress noisy access logs for static files; show proxy calls
        if '/horizons-proxy' in (args[0] if args else ''):
            print(f'[Horizons proxy] {args[0]}')

if __name__ == '__main__':
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    print(f'Artemis II proxy server running at http://localhost:{PORT}')
    print(f'Proxying JPL Horizons via /horizons-proxy')
    http.server.HTTPServer(('', PORT), Handler).serve_forever()
