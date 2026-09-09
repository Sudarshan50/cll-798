import React from "react";

export function ErrorBanner() {
  return (
    <div className="err">
      <b>The atlas could not load.</b>
      <span>Reload the page to try again.</span>
      <span className="hint">
        It needs a current desktop browser with WebGL2 and a stable connection —
        the compound layer is a 292 MB download.
      </span>
    </div>
  );
}

export function ModeToggle({ mode, onChange }) {
  return (
    <div className="modes">
      <button className={mode === "network" ? "on" : ""}
        onClick={() => onChange("network")}>Network</button>
      <button className={mode === "compounds" ? "on" : ""}
        onClick={() => onChange("compounds")}>Network + compounds</button>
    </div>
  );
}

/** Bottom-right stack. Hint and toggle used to be two independently positioned
 *  fixed elements whose offsets were hand-tuned, so they overlapped. One flex
 *  column keeps them apart whatever either one contains. */
export function CornerControls({ labels, onChange }) {
  return (
    <div className="corner">
      <div className="hint">
        <span><b>left click</b> highlight</span>
        <span><b>right click</b> details</span>
        <span><b>drag</b> pan · <b>scroll</b> zoom</span>
      </div>
      <div className="ctl">
        <label><input type="checkbox" checked={labels}
          onChange={(e) => onChange(e.target.checked)} /> labels</label>
      </div>
    </div>
  );
}

export function SelectionBar({ sel, onFit, onClear }) {
  return (
    <div className="selbar">
      <span className="dot" style={{ background: sel.col }} />
      <div className="who"><div className="n">{sel.n}</div></div>
      <button className="fit" onClick={onFit}
        title="Frame this class and everything linked to it">fit</button>
      <button className="x" onClick={onClear} aria-label="Clear selection">×</button>
    </div>
  );
}

export function Unsupported({ minW, minH }) {
  return (
    <div className="gate">
      <div className="gate-in">
        <h1>Desktop only</h1>
        <p>
          The atlas draws 4,824 taxonomy classes and streams 292&nbsp;MB of compound
          geometry into a WebGL2 buffer. A phone or tablet cannot allocate that, and
          the class details open on right-click, which touch has no equivalent for.
        </p>
        <p className="req">
          Needs a mouse or trackpad and a window of at least
          {" "}<b>{minW}&nbsp;&times;&nbsp;{minH}&nbsp;px</b>.
        </p>
        <p className="alt">
          Open this page on a laptop or desktop to explore the atlas.
        </p>
      </div>
    </div>
  );
}
