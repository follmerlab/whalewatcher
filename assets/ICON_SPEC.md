# whalewatcher application icon

**Illustration type:** `composite` (flat vector mark, not a publication figure)
**Canvas:** 512x512, `viewBox="0 0 512 512"`
**Context:** application icon — Tk window icon, favicon, GitHub social preview

## Concept

The repo name is a pun: ORCA is a killer whale, so a tool that reads ORCA output is a
whale watcher. The mark takes that literally — **the badge is the view through a
telescope**. The circle is the eyepiece field, framed by the field stop, with edge
vignette and a glass glint. No reticle or crosshair.

The sea is drawn as a sine profile, a quiet nod to the vibrational-mode tab.

Self-contained circular badge, so it needs no help from the host background and works
unchanged on white and on the app's `#1a1a2e` canvas.

## Why the scope framing replaced the binoculars

Two earlier passes drew the watcher as an object inside the scene — first a pair of eyes,
then binoculars. Making the frame itself the instrument turned out to be strictly better
on three counts:

1. **It frees the field.** With no binoculars competing for the top third, the whale is
   drawn ~15% larger on screen despite the field being tighter.
2. **The barrel ring reads at 16 px.** A dark rim gives the icon a crisp boundary where a
   flat cyan disc just dissolves into the page.
3. **It removed the need for two artworks.** The binoculars version needed a separate
   reduced artwork in the 16/32 slots; the scope version carries all six sizes itself.

## Legibility, measured not asserted

Judged from nearest-neighbour blow-ups of the real rasters, on both backgrounds.

| Size | `scope` (shipped) | `small` | `detailed` (binoculars) |
|---|---|---|---|
| 128 px | excellent | excellent | excellent |
| 64 px | excellent | very clean | good |
| 48 px | excellent | good | good |
| 32 px | **good** — fin, belly, eye patch all read | good | busy; whale is a smudge |
| 16 px | **best of the three** — dark rim, light field, dark body mass | pale disc, weak boundary | fails, indistinct blur |

**Honest note on 16 px:** nobody will identify a killer whale at 16 px in any version.
`scope` degrades to a dark-rimmed lens with a dark shape in it — distinctive and crisp
rather than mushy, which was the bar. The rim is what carries it.

`whalewatcher.ico` uses `scope` in all six slots. Verified to contain six frames.

## Palette

| Role | Hex |
|---|---|
| Field / sky | `#8ED8EC` |
| Orca body | `#101A28` |
| Belly, eye patch | `#FFFFFF` |
| Field stop / barrel | `#22384E` → `#111F2E` → `#060D16` |
| Vignette | `#04101C`, 0 → 0.86 alpha |
| Sea, near / far | `#2A8CB0` / `#3FA8CC` |
| Catchlight | `#7FB4CE` at 38% |

Gradients are used only for the barrel, vignette, and glint. All three are large-scale
luminance ramps, which downsample cleanly; the rule being avoided is gradients carrying
*fine detail*, which is what muds up at small sizes.

## Design decisions and what was rejected

Each of these was rendered and looked at before being cut.

- **Eye patch as a two-lobed p orbital** (the original brief's preferred fusion) — cut.
  Two white ovals side by side on a head read as googly eyes, not an orbital. Replaced
  with a single anatomically-placed eye patch.
- **Symmetric watching eyes** — cut from the primary. A symmetric pair inside a circle
  combines with the badge into an unmistakable smiley face, with the orca as the mouth.
  Survives in `whalewatcher_icon_eyes.svg` only because the eyes were moved off-centre
  and rotated, which breaks the gestalt.
- **Grey saddle patch** — cut. At every size it read as a chip out of the whale's back.
- **Vertical forked tail** — cut. That is a fish caudal fin; whale flukes are horizontal.
  Reworked into a broad flat fluke, which is most of why the animal stopped reading as a
  generic fish.
- **Flank blaze** — cut. Merged with the belly into one white mass running into the
  flukes, making the tail look detached.
- **Traced belly outline** — cut. Its upper edge sat below the silhouette's lower edge, so
  the clip ate nearly all of it. Replaced with an oversized ellipse that cuts across the
  body at y~295.
- **Thick lit metal ring** — cut. The first scope pass read as a ship's porthole, which is
  the wrong instrument. A telescope eyepiece is a narrow, dark, nearly flat aperture with
  heavy edge falloff, so the ring was thinned and flattened and the vignette deepened.
- **Butt-capped catchlight arc** — cut. It terminated abruptly and left a visible stub
  notch on the right edge. Round caps, both ends buried in the vignette.
- **Deep mid-range vignette** — cut. It swallowed the tail flukes. The falloff now stays
  clear to 78% of the radius and then ramps hard, keeping the tube read without eating
  the subject.
- **Crosshairs / reticle** — never added. Explicitly excluded by the brief, and it would
  have read as a gunsight, which is a bad look for a whale.

## Variants

| File | Status |
|---|---|
| `whalewatcher_icon.svg` | **canonical master and shipped artwork** — telescope framing, used for every `.ico` size |
| `variants/whalewatcher_icon_binoculars.svg` | superseded; watcher drawn as an object in the scene |
| `variants/whalewatcher_icon_eyes.svg` | superseded; playful rather than professional |
| `variants/whalewatcher_icon_small.svg` | flat mark, no scope framing; kept for contexts wanting one |

Variants are not built by default. `make_icon.py binoculars` (or `eyes` / `flat`) renders one
on demand.

## Regenerating

**Not needed to run the app.** Every raster is committed — `whalewatcher.ico` and
`whalewatcher_icon_{16,32,48,64,128,256,512}.png` — and `orca_vib_viewer.py` loads those
files directly. Clone and run; nothing here is a build step.

`make_icon.py` is a maintainer tool, only for regenerating the rasters after editing
`whalewatcher_icon.svg`. It is the only thing that needs Playwright and Chromium, which are
deliberately *not* project dependencies.

```bash
python make_icon.py
```

Renders each master once at 1024 px through headless Chromium, then Lanczos downsamples.
Downsampling from a high-res render beats rasterising the vector natively at 16 px, which
throws away the antialiasing that makes a small icon readable.

`cairosvg` was specified in the brief but cannot load libcairo on this machine, so
Playwright does the rasterising and Pillow handles the `.ico` and contact sheet. Both are
already dependencies of the project's tooling; the pipeline is fully reproducible.

## Wiring into the app

`OrcaVibViewer._set_window_icon()` applies it, called from `__init__`:

1. **All platforms** — `iconphoto(False, *pngs)` with the 256/128/64/48/32/16 PNGs. The
   `default` flag is the load-bearing detail: `iconphoto(True, ...)` and
   `iconbitmap(default=...)` set only the application/class default and leave the window's
   own icon slots empty (`WM_GETICON` returns null). `default=False` issues `WM_SETICON`
   against this window. Measured directly:

   | call | `WM_GETICON` |
   |---|---|
   | `iconphoto(True, ...)` | `(0, 0, 0)` |
   | `iconphoto(False, ...)` | real handles, distinct for small and big |

2. **Windows, additionally** — `iconbitmap(default=...)` with the `.ico`, applied *after*
   step 1 (verified not to disturb it), so message boxes and future toplevels inherit it.
3. **Assets missing** — every failure path falls through silently to Tk's default feather.
   A missing icon must never stop the viewer opening.

`set_windows_app_id()` runs before `super().__init__()`, since the taskbar reads process
identity at window-creation time. It gives the app its own taskbar button instead of
grouping under `python.exe`. It does **not** fix the taskbar icon — see below.

One non-obvious detail: Tk does not retain a reference to `PhotoImage` objects. They are
held on `self._icon_images`; without that they are garbage collected and the icon silently
reverts to the default.

### Known limitation: Windows taskbar shows the Tk feather

Title bar and Alt-Tab are correct. The taskbar button is not, and it is not for want of
trying. What was established by measurement:

- the window's `WM_GETICON` handles were extracted and rendered — they contain the orca
- the window class `GCLP_HICON` / `GCLP_HICONSM` likewise contain the orca
- enumerating the process's windows found exactly one taskbar-eligible window, carrying
  those correct icons — no stray Tk window owning the button
- `python.exe`'s own embedded icon is the Python logo, not a feather, so the taskbar is not
  simply showing the interpreter's icon

Things tried that changed nothing: raw win32 `LoadImageW` + `WM_SETICON`, overwriting the
class icon with `SetClassLongPtr`, clearing the AppUserModelID, and setting a fresh unused
one.

That exhausts what a process can do to its own windows, which points at the shell icon
cache — plausibly keyed on the interpreter path, so every Tk script on this interpreter
inherits whatever was cached first. The durable fixes both live outside the process:
launch via a `.lnk` carrying the icon and a matching AppUserModelID, or freeze to an `.exe`
with the icon as a resource. Neither is worth it for a lab tool right now. Tracked as
[issue #17](https://github.com/follmerlab/whalewatcher/issues/17).

The `iconphoto(False, ...)` call is still doing real work regardless — it is what makes the
icon correct on Linux and macOS, where `.ico` is not supported at all.

Earlier paths were exercised against a live Tk root — the `.ico` branch, the PNG
fallback (forced by simulating a Tk that rejects `.ico`, confirming all six images load),
and the missing-assets branch.
