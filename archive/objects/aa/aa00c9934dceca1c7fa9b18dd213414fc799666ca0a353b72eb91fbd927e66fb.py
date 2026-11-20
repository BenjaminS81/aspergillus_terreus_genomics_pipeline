def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda: f.read(1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def write_tiff(page, path):
    """Baseline RGB Deflate TIFF, strip by strip; never allocate a full page."""
    raster_bounds = (page.rect * pymupdf.Matrix(S, S)).irect
    width = raster_bounds.width
    height = raster_bounds.height
    display = page.get_displaylist()
    offsets = []
    counts = []
    strip_hashes = {}
    pending = path.with_suffix('.pending.tif')
    with pending.open('w+b') as h:
        h.write(b'II' + struct.pack('<HI', 42, 0))
        for i, y in enumerate(range(0, height, ROWS)):
            bottom = min(y + ROWS, height)
            clip = pymupdf.Rect(0, max(0, (y - 1) / S), page.rect.width, min(page.rect.height, (bottom + 1) / S))
            pix = display.get_pixmap(matrix=pymupdf.Matrix(S, S), colorspace=pymupdf.csRGB, alpha=False, clip=clip)
            pixels = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, 3)
            band = np.ascontiguousarray(pixels[y - pix.y:bottom - pix.y, :width, :])
            assert band.shape == (bottom - y, width, 3), (i, band.shape, width, height, pix.irect)
            raw = band.tobytes()
            compressed = zlib.compress(raw, 4)
            offsets.append(h.tell())
            counts.append(len(compressed))
            h.write(compressed)
            if i in [0, height // ROWS // 2, (height - 1) // ROWS]:
                strip_hashes[i] = hashlib.sha256(raw).hexdigest()
            del raw, compressed, pix, pixels, band

        def payload(data):
            if h.tell() % 2:
                h.write(b'\x00')
            off = h.tell()
            h.write(data)
            return off
        bits = payload(struct.pack('<HHH', 8, 8, 8))
        strip_offsets = payload(struct.pack('<' + 'I' * len(offsets), *offsets))
        strip_counts = payload(struct.pack('<' + 'I' * len(counts), *counts))
        res = payload(struct.pack('<II', DPI, 1))
        software = b'Terreus print export v2; rendered from vector PDF at 600 dpi\x00'
        soft = payload(software)
        if h.tell() % 2:
            h.write(b'\x00')
        ifd = h.tell()
        tags = [(256, 4, 1, width), (257, 4, 1, height), (258, 3, 3, bits), (259, 3, 1, 8), (262, 3, 1, 2), (273, 4, len(offsets), strip_offsets), (277, 3, 1, 3), (278, 4, 1, ROWS), (279, 4, len(counts), strip_counts), (282, 5, 1, res), (283, 5, 1, res), (284, 3, 1, 1), (296, 3, 1, 2), (305, 2, len(software), soft)]
        h.write(struct.pack('<H', len(tags)))
        for tag, typ, n, value in tags:
            h.write(struct.pack('<HHII', tag, typ, n, value))
        h.write(struct.pack('<I', 0))
        h.seek(4)
        h.write(struct.pack('<I', ifd))
    Image.MAX_IMAGE_PIXELS = None
    with Image.open(pending) as im:
        assert im.format == 'TIFF' and im.size == (width, height) and (im.mode == 'RGB')
        assert im.info['dpi'] == (600.0, 600.0), im.info
        assert im.tag_v2[259] == 8
    with pending.open('rb') as h:
        for i, digest in strip_hashes.items():
            h.seek(offsets[i])
            raw = zlib.decompress(h.read(counts[i]))
            assert hashlib.sha256(raw).hexdigest() == digest
    os.replace(pending, path)
    return {'tiff': 'print_600dpi/' + path.name, 'dpi': DPI, 'pixels': [width, height], 'tiff_sha256': sha(path), 'tiff_bytes': path.stat().st_size, 'sampled_strips_verified': True}
