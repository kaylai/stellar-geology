# Contributing a Gallery Example

Gallery examples are Python scripts in the `examples/` directory. At build time, [sphinx-gallery](https://sphinx-gallery.github.io) runs each script, captures its output and figures, and renders it as a docs page with downloadable `.py` and `.ipynb` versions. Everything under `docs/auto_examples/` is generated — never edit it directly.

To add an example:

1. Create `examples/plot_<your_example>.py`. The `plot_` prefix is required — scripts without it are not executed.
2. Start the file with a docstring containing the title, underlined with `=`:

   ```python
   """
   My Example Title
   ================

   One or two sentences describing the example.
   """
   ```

3. Split the script into cells with `# %%`. Comment lines directly after a `# %%` are rendered as prose (reStructuredText) between the code blocks:

   ```python
   # %%
   # Build a star and compute its composition.

   import stellar_geology as sg

   star = sg.Star(stellar_dex={"Fe": 0.1})
   ```

4. Build the docs and check your example renders and runs cleanly:

   ```bash
   cd docs
   make html
   ```

That's it — the example appears in the gallery automatically. See [plot_unit_conversions.py](https://github.com/kaylai/stellar-geology/blob/main/examples/plot_unit_conversions.py) for a complete reference.
