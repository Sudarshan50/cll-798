import React from "react";

export function ErrorBanner({ message }) {
  return <div className="err"><b>Could not start.</b><span>{message}</span></div>;
}

export function ModeToggle({ mode, onChange }) {
  return (
    <div className="modes">
      <button className={mode === "network" ? "on" : ""}
        onClick={() => onChange("network")}>1 · Network</button>
      <button className={mode === "compounds" ? "on" : ""}
        onClick={() => onChange("compounds")}>3 · Network + compounds</button>
    </div>
  );
}

export function LabelToggle({ labels, onChange }) {
  return (
    <div className="ctl">
      <label><input type="checkbox" checked={labels}
        onChange={(e) => onChange(e.target.checked)} /> labels</label>
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

export function Hint() {
  return (
    <div className="hint">
      <span><b>left click</b> highlight</span>
      <span><b>right click</b> details</span>
      <span><b>drag</b> pan · <b>scroll</b> zoom</span>
    </div>
  );
}
