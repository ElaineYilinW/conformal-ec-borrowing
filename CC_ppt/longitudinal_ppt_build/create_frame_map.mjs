import fs from "node:fs";
import path from "node:path";

const out = path.resolve("CC_ppt/longitudinal_ppt_build/template-frame-map.json");
const changed = {
  25: [
    ["3", "rewrite"], ["5", "rewrite-and-reposition"], ["6", "rewrite-and-reposition"],
    ["18", "replace"], ["8", "rewrite-and-reposition"], ["9", "rewrite-and-reposition"],
    ["19", "replace"], ["11", "rewrite-and-reposition"], ["20", "replace"],
    ["13", "rewrite-and-reposition"], ["14", "rewrite-and-reposition"],
    ["21", "replace"], ["16", "rewrite-and-reposition"], ["17", "keep"],
  ],
  26: [
    ["3", "rewrite"], ["9", "replace"], ["10", "replace"],
    ["7", "rewrite-and-reposition"], ["8", "keep"],
  ],
  27: [
    ["3", "rewrite"], ["5", "rewrite-and-reposition"],
    ["8", "replace"], ["7", "keep"],
  ],
  28: [
    ["3", "rewrite"], ["5", "rewrite-and-reposition"],
    ["9", "replace"], ["7", "rewrite-and-reposition"], ["8", "keep"],
  ],
};

const roles = {
  25: "simulation DGP and analysis model",
  26: "simulation evidence: two-sided screen in one representative draw",
  27: "simulation comparison: whole-trajectory score detection",
  28: "simulation summary table and operating characteristics",
};

const outputSlides = Array.from({ length: 34 }, (_, index) => {
  const slide = index + 1;
  return {
    outputSlide: slide,
    sourceSlide: slide,
    narrativeRole: roles[slide] || "preserve source slide unchanged",
    reuseMode: "duplicate-slide",
    editTargets: (changed[slide] || []).map(([sourceElementId, action]) => ({
      sourceElementId,
      action,
    })),
  };
});

fs.mkdirSync(path.dirname(out), { recursive: true });
fs.writeFileSync(out, `${JSON.stringify({ outputSlides }, null, 2)}\n`);
console.log(out);

