# 🎈 Blank app template

A simple Streamlit app template for you to modify!

[![Open in Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://blank-app-template.streamlit.app/)

### How to run it on your own machine

Prerequisite: install `uv` if you don't already have it.

```
$ curl -LsSf https://astral.sh/uv/install.sh | sh
```

1. Sync the dependencies

   ```
   $ uv sync
   ```

2. Run the app

   ```
   $ uv run streamlit run streamlit_app.py
   ```

### Primer design

Add sequences in assembly order. The app designs primers at every adjacent
sequence junction by searching for A-to-T overhangs 6–12 bases long and
preferring the selected overhang length. Each junction displays the selected
overhang, sequence ends, and forward/reverse primers in 5′→3′ orientation.

Set the primer binding length and preferred overhang length in the sidebar.
The overlap preference on the left-hand sequence controls whether a junction
searches that sequence's end, the next sequence's start, or both.
