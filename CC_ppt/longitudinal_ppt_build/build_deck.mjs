import fs from "node:fs/promises";
import path from "node:path";
import { FileBlob, PresentationFile } from "@oai/artifact-tool";

const root = "/Users/yilinwang/Desktop/Herb";
const starterPptxPath = process.env.STARTER_PPTX || path.join(process.cwd(), "template-starter.pptx");
const outputPptxPath = process.env.OUTPUT_PPTX || path.join(root, "CC_ppt", "CC_presentation_longitudinal_sim.pptx");

const deck = await PresentationFile.importPptx(await FileBlob.load(starterPptxPath));

const COLORS = {
  blue: "#2F7DC8",
  navy: "#1F2A44",
  body: "#333740",
};

const pos = ([left, top, width, height]) => ({ left, top, width, height });
const shapeByName = (slide, name) => {
  const found = slide.shapes.items.find((item) => item.name === name);
  if (!found) throw new Error(`Missing shape ${name}`);
  return found;
};
const imageByName = (slide, name) => {
  const found = slide.images.items.find((item) => item.name === name);
  if (!found) throw new Error(`Missing image ${name}`);
  return found;
};
const setPlain = (shape, text, frame, style = {}) => {
  shape.text = text;
  shape.position = pos(frame);
  shape.text.style = {
    fontSize: 18,
    typeface: "Arial",
    color: COLORS.body,
    autoFit: "shrinkText",
    verticalAlignment: "middle",
    insets: { left: 0, right: 0, top: 0, bottom: 0 },
    ...style,
  };
};
const setRich = (shape, paragraphs, frame, fontSize = 18) => {
  shape.text = paragraphs;
  shape.position = pos(frame);
  shape.text.style = {
    fontSize,
    typeface: "Arial",
    color: COLORS.body,
    autoFit: "shrinkText",
    verticalAlignment: "top",
    insets: { left: 0, right: 0, top: 0, bottom: 0 },
  };
};
const replaceImage = async (image, filePath, frame, alt) => {
  const bytes = await fs.readFile(filePath);
  image.replace({ blob: bytes, contentType: "image/png", alt, fit: "contain" });
  image.frame = pos(frame);
  image.crop = { left: 0, top: 0, right: 0, bottom: 0 };
  image.lockAspectRatio = false;
};
const replaceTitle = (slide, oldText, newText) => {
  const title = shapeByName(slide, "TextBox 2");
  title.text.get(oldText).text = newText;
};
const para = (runs, opts = {}) => ({ runs, spaceAfter: opts.spaceAfter ?? 8 });
const run = (text, opts = {}) => ({
  run: text,
  textStyle: {
    bold: opts.bold,
    color: opts.color || COLORS.body,
    fontSize: `${opts.fontSize || 18}pt`,
    typeface: "Arial",
  },
});

// Slide 25 — full DGP.
{
  const slide = deck.slides.getItem(24);
  replaceTitle(slide, "Simulation: the data-generating process", "12-visit longitudinal DGP: shared covariate effects");

  setPlain(shapeByName(slide, "TextBox 4"),
    "Samples and mild covariate shift (S=1 trial; S=0 external; centered binary terms use B tilde):",
    [52.8, 142, 1180, 32], { fontSize: 18.5 });
  shapeByName(slide, "Rectangle 5").position = pos([55.7, 174, 1167, 74]);
  await replaceImage(imageByName(slide, "Picture 6"),
    path.join(root, "CC_ppt", "eq", "eq_long_sample_source.png"),
    [67, 181, 1143, 60], "Sample sizes and logistic source model");

  setPlain(shapeByName(slide, "TextBox 7"),
    "Baseline construction, then one shared covariate coefficient vector and a small linear time slope:",
    [52.8, 253, 1180, 32], { fontSize: 18.5 });
  shapeByName(slide, "Rectangle 8").position = pos([55.7, 285, 1167, 66]);
  await replaceImage(imageByName(slide, "Picture 9"),
    path.join(root, "CC_ppt", "eq", "eq_long_baseline.png"),
    [67, 294, 1143, 48], "Baseline covariate and baseline outcome distribution");

  shapeByName(slide, "Rectangle 10").position = pos([55.7, 358, 1167, 106]);
  await replaceImage(imageByName(slide, "Picture 11"),
    path.join(root, "CC_ppt", "eq", "eq_long_outcome.png"),
    [67, 366, 1143, 90], "Twelve-visit outcome mean and compound-symmetric covariance");

  setPlain(shapeByName(slide, "TextBox 12"),
    "Thirty percent of external controls are incompatible; the discrepancy has a covariate-dependent level plus a small time drift:",
    [52.8, 473, 1180, 32], { fontSize: 18.5 });
  shapeByName(slide, "Rectangle 13").position = pos([55.7, 505, 1167, 73]);
  await replaceImage(imageByName(slide, "Picture 14"),
    path.join(root, "CC_ppt", "eq", "eq_long_contam.png"),
    [67, 513, 1143, 56], "External-control contamination model");

  setPlain(shapeByName(slide, "TextBox 15"),
    "Analysis: nlme::gls(Y ~ Y0 + X1 + X2 + B1 + ... + B4 + time, corCompSymm(~1|id), REML). All 12 visits estimate the shared coefficients used in visit-12 AIPW. Estimand: trial-population ATE at visit 12; true tau = 0.70. Residual variance = 4 and CS correlation = 0.30.",
    [52.8, 585, 1098, 90], { fontSize: 17.2, bold: true, color: COLORS.navy, verticalAlignment: "top" });
}

// Slide 26 — representative draw and explicit two-sided weighted p-value.
{
  const slide = deck.slides.getItem(25);
  replaceTitle(slide, "The screen in action (one draw at the reference cell)", "The two-sided screen in action (one representative draw)");

  const body = [
    para([run("Reference draw. ", { bold: true, color: COLORS.blue, fontSize: 19 }), run("Seed 2 gives 34 contaminated ECs; the primary screen detects 94.1% and wrongly excludes 9.1% of clean ECs.")], { spaceAfter: 10 }),
    para([run("Top. ", { bold: true, color: COLORS.blue, fontSize: 19 }), run("After the RCT-control 12-visit GLS, compatible EC residuals track zero, while incompatible trajectories sit roughly 6–8 units lower across visits.")], { spaceAfter: 10 }),
    para([run("Bottom. ", { bold: true, color: COLORS.blue, fontSize: 19 }), run("Bad-EC p-values concentrate below gamma=0.10; clean-EC p-values remain spread over [0,1].")], { spaceAfter: 10 }),
    para([run("Cross-fitting. ", { bold: true, color: COLORS.blue, fontSize: 19 }), run("Five folds refit both the shared-coefficient GLS and the logistic source model; weighted RCT-control ranks are compared with each EC.")], { spaceAfter: 10 }),
    para([run("Two-sided by construction. ", { bold: true, color: "#C74440", fontSize: 19 }), run("The nonconformity score uses an absolute residual contrast, so positive and negative incompatibility are treated symmetrically.")]),
  ];
  setRich(shapeByName(slide, "TextBox 6"), body, [52.8, 145, 566, 335], 18);

  await replaceImage(imageByName(slide, "Picture 4"),
    path.join(root, "CC_ppt", "longitudinal_assets", "one_draw_two_sided.png"),
    [650, 145, 575, 516], "Representative residual trajectories and two-sided p-value histogram");
  await replaceImage(imageByName(slide, "Picture 5"),
    path.join(root, "CC_ppt", "eq", "eq_long_score_p.png"),
    [55, 492, 565, 155], "Cross-fitted average residual score and weighted two-sided conformal p-value");
}

// Slide 27 — actual 200-replication score comparison.
{
  const slide = deck.slides.getItem(26);
  replaceTitle(slide, "When can the screen detect contamination?", "Whole-trajectory scores preserve the detectable signal");
  const body = [
    para([run("Why average works. ", { bold: true, color: COLORS.blue, fontSize: 19 }), run("The incompatible DGP is mainly a common level discrepancy plus 0.10 × time. Averaging all visits keeps that level signal; the visit-12-minus-visit-1 contrast cancels it and leaves only 1.1 units.")], { spaceAfter: 9 }),
    para([run("Signal-to-noise. ", { bold: true, color: COLORS.blue, fontSize: 19 }), run("Standardized SNR is 5.37 for the visit average, 3.49 for visit 12 alone, and only 0.465 for the first-to-last contrast.")], { spaceAfter: 9 }),
    para([run("Observed screen. ", { bold: true, color: COLORS.blue, fontSize: 19 }), run("At the two-sided gamma=0.10 cutoff, average/internal-CS detects 93.4% with 8.3% wrong exclusion; final-only detects 85.2%, but c1 detects just 11.7%.")], { spaceAfter: 9 }),
    para([run("Covariance choice. ", { bold: true, color: COLORS.blue, fontSize: 19 }), run("Estimate score covariance from internal RCT controls. Estimating it from contaminated ECs inflates the scale and masks the signal.")]),
  ];
  setRich(shapeByName(slide, "TextBox 4"), body, [52.8, 145, 1171, 208], 18);
  await replaceImage(imageByName(slide, "Picture 5"),
    path.join(root, "CC_ppt", "longitudinal_assets", "score_comparison.png"),
    [80, 360, 1120, 300], "Detection and wrong-exclusion rates for longitudinal score variants");
}

// Slide 28 — main 200 x 500 operating characteristics.
{
  const slide = deck.slides.getItem(27);
  replaceTitle(slide, "Result: beats RCT-only on all four criteria at once", "Result: essentially unbiased, lower RMSE, correct type-I, higher power");
  shapeByName(slide, "Rectangle 4").position = pos([115, 140, 1050, 390]);
  await replaceImage(imageByName(slide, "Picture 5"),
    path.join(root, "CC_ppt", "eq", "tab_long_main.png"),
    [130, 150, 1020, 370], "Main longitudinal simulation operating characteristics");

  const body = [
    para([run("Primary conformal average. ", { bold: true, color: COLORS.blue, fontSize: 19 }), run("Bias = 0.016; RMSE = 0.285, an 11.9% reduction from RCT-longitudinal (0.323); type-I = 0.055; power rises from 0.555 to 0.715.")], { spaceAfter: 8 }),
    para([run("Screen quality. ", { bold: true, color: COLORS.blue, fontSize: 19 }), run("Detection = 0.934 and wrong exclusion = 0.083; on average 64.1 clean and 2.0 bad ECs remain, with borrowing coefficient 0.507.")], { spaceAfter: 8 }),
    para([run("Monte Carlo inference. ", { bold: true, color: COLORS.blue, fontSize: 19 }), run("200 replications and 500 full-pipeline subject bootstraps. Type-I and power use two-sided 95% percentile intervals; the conformal p-value is also two-sided through |V|.")]),
  ];
  setRich(shapeByName(slide, "TextBox 6"), body, [52.8, 536, 1100, 132], 17.4);
}

const pptx = await PresentationFile.exportPptx(deck);
await pptx.save(outputPptxPath);
console.log(outputPptxPath);

