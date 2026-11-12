import React from "react";

export function MethodTab() {
  return (
    <div className="x-prose">
      <h3>Phase 2 · a confusability overlay</h3>
      <p>
        ClassyFire gives every compound one direct parent and, on average, 12.6 alternative
        parents. Two classes that keep being assigned to the same compounds are linked; the
        link is weighted by pointwise mutual information, so frequent classes do not win
        just by being large. Only <em>competing</em> pairs are kept — those that look like
        two answers to the same question rather than one class describing another.
      </p>
      <p>
        Each class gets an ambiguity score <code>H</code>: the entropy of its competing
        partners, filtered (overlay only) and unfiltered (all alternative parents). The
        overlay is then tested three ways: are the classes the paper's experts marked as
        wrong more ambiguous than size-matched controls (H1a); is the graph more clustered
        and modular than a degree-preserving rewiring (H1b); and do its communities recut
        the SuperClass partition rather than redraw it (H1c).
      </p>

      <h3>Bipartite graph and projections</h3>
      <p>
        The labels form a bipartite graph: compounds on one side, classes on the other, an
        edge for every label. Projecting it onto the classes links two classes that share
        a compound; projecting onto compounds links two compounds that share a class. The
        backbone drawn on the map keeps only class–class links that occur more often than
        the bipartite configuration model (BiCM) predicts from the degrees alone.
      </p>

      <h3>Caveats</h3>
      <ul>
        <li><b>Ambiguity is not error.</b> A high <code>H</code> says a class sits among
          competing alternatives, not that its assignments are wrong.</li>
        <li><b>n = 20.</b> The only ground truth is 20 expert-flagged classes (from 30 FP/FN
          markings). That is a pilot; one class moves the p-value.</li>
        <li><b>Alternative parents are descriptors by design.</b> Most co-assignments are
          intended, so the descriptive-vs-competing filter carries the result and is swept
          for robustness.</li>
        <li><b>Projection inflates clustering.</b> Every class becomes a clique in the
          compound projection and every compound a clique among its classes, so raw
          clustering is high by construction. Compare against a null, never in absolute terms.</li>
        <li><b>BiCM is the null.</b> It keeps every node's degree on both sides and asks
          what co-occurrence is left; edges that survive it are the ones worth reading.</li>
      </ul>

      <h3>Where the numbers come from</h3>
      <p>
        Full write-ups: <code>docs/phase2_results.md</code> and{" "}
        <code>docs/bipartite_extension.md</code>; design and pre-registration in{" "}
        <code>docs/phase2_plan.md</code>. The panel reads <code>phase2.json</code> and{" "}
        <code>bipartite.json</code>; anything they leave out is shown as “not computed”.
      </p>
    </div>
  );
}
