"""Saved-input report figure: reference_context."""
from ._common import Path, frame, save

def plot(source, output, resources, image_path=None):
    import shutil
    import matplotlib.image as mpimg
    image = Path(image_path) if image_path else resources / 'reference/reference_context.png'
    out = Path(output)
    out.mkdir(parents=True, exist_ok=True)
    pixels = mpimg.imread(image)
    fig, axes = frame('7NFB deposited asymmetric-unit context — reference 7NFB|B|1', size=(12, 10))
    axes[0].imshow(pixels)
    axes[0].axis('off')
    save(fig, out, '05_reference_context', 'Full original 7NFB image and caption, not a substituted selected complex: selected receptor blue, other polymers grey, peptides orange and GEN green across receptor copies. No LaTeX crop applied. Offline route reuses this genuine saved render; render_reference verifies an actual new PyMOL execution separately.')
    shutil.copy2(image, out / '05_reference_context.png')
