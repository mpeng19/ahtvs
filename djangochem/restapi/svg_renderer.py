from rest_framework import renderers
import svgwrite

from molgen import blockdraw
from blocks import block
from django.conf import settings
import os

IMG_CACHE = getattr(settings, 'SVG_CACHE_DIR', '~/mol_img_cache')
IMG_CACHE = os.path.expanduser(IMG_CACHE)

DEFAULT_SVG_SIZE = getattr(settings, 'SVG_SIZE', 200)
DEFAULT_SVG_COLS = getattr(settings, 'SVG_COLS', 5)
GROUPS = getattr(settings, 'SVG_ATOM_GROUPS', [])  # examples: "SO3H,COOH,PO3H2" or "CF3"


class SVGRenderer(renderers.BaseRenderer):
    media_type = 'image/svg+xml'
    format = 'svg'
    charset = 'utf8'
    render_style = 'binary'

    def _svg_text_from_data(self, data, refresh=False):
        filepath = blockdraw.write_svg(block.Block(data["smiles"]),
                                       path=IMG_CACHE,
                                       force=refresh,
                                       groups=GROUPS)
        return open(filepath, 'r').read()

    def render(self, data, media_type=None, renderer_context=None):
        request = renderer_context['request']
        try:
            size = int(request.query_params.get('svgsize', DEFAULT_SVG_SIZE))
        except ValueError:
            size = DEFAULT_SVG_SIZE

        try:
            cols = int(request.query_params.get('svgcols', DEFAULT_SVG_COLS))
        except ValueError:
            cols = DEFAULT_SVG_COLS

        refresh = request.query_params.get('svgrefresh', False)
        if refresh in ['True', 'true', '1']:
            refresh = True
        else:
            refresh = False
        if "results" in data:
            results = data["results"]
            count = len(results)
            if count == 1:
                return self._svg_text_from_data(results[0], refresh)
            max_rows = int(count / cols) + 1
            doc_size = ("{}px".format(size*cols), "{}px".format(size*max_rows))
            svg_document = svgwrite.Drawing(size=doc_size)
            for i, result in enumerate(results):
                row = int(i / cols)
                col = i % cols
                href = request.path_info + str(result["id"])+"?format=svg"
                if refresh:
                    href += "&svgrefresh=true"
                svg_document.add(svg_document.image(href=href, insert=(col*size, row*size), size=(size, size)))

            return svg_document.tostring()
        else:
            return self._svg_text_from_data(data, refresh)
