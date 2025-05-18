import streamlit as st


def inject_base_css() -> None:
    """Put here every <style> block that should be loaded once."""
    st.markdown(
        """
        <style>
            /* logo + page width tweaks, tooltip css, etc. */
            .tooltip      {display:inline-block; position:relative; font-size:15px;}
            .tooltiptext  {
                visibility:hidden; width:340px; background:#f9f9f9; color:#333;
                text-align:left; padding:6px 12px; border:2px solid #ccc;
                border-radius:6px; position:absolute; z-index:1;
                top:125%; left:50%; transform: translateX(-50%);
                box-shadow:0 2px 8px rgba(0,0,0,0.15); opacity:0; transition:opacity .25s; top:-1px; 
            }
            .tooltip:hover .tooltiptext {visibility:visible; opacity:1;}

            div.row-widget.stButton {margin:0rem;}
            div[data-testid="stHorizontalBlock"] {margin-bottom:0.5rem!important;}
            .block-container {max-width:1100px;}
        </style>
        """,
        unsafe_allow_html=True,
    )

box_style = """
<div style="
    border: 2px solid {border};
    border-radius: 10px;
    padding: 20px;
    margin: 10px;
    background-color: {bg};
    text-align: center;">
    <h3 style="color:{title_color};">{title}</h3>
    <p style="color:#000000;">{content}</p>
</div>
"""


def html_box(title: str, content: str,
             border: str = "#28a745",
             bg: str = "#f9f9f9",
             fg: str = "#000000") -> str:
    """Return a ready-to-embed HTML info box."""
    title_color = "#000000" if bg.lower() != "#f9f9f9" else border
    return box_style.format(
            title=title,
            content=content,
            border=border,
            bg=bg,
            title_color=title_color
    )

def tooltip_icon(text: str, body_html: str) -> str:
    """Embed the little ℹ️ with a custom tooltip bubble."""
    return f"""
    <span class="tooltip">{text}
        <span class="tooltiptext">{body_html}</span>
    </span>
    """

def show_info_box(title: str, content: str, **kwargs) -> None:
    # allow color=… as an alias for bg=…
    if "color" in kwargs:
        kwargs["bg"] = kwargs.pop("color")

    border = kwargs.get("border", "#28a745")
    bg     = kwargs.get("bg", "#f9f9f9")
    #fg     = kwargs.get("fg", "#000000")

    st.markdown(html_box(title, content, border, bg),
                unsafe_allow_html=True)


MISSING_INPUT_CSS = """
<style>
.tooltip-wrapper {position:relative; display:inline-block; width:100%;}
.tooltip-bubble  {
    display:none; position:absolute; top:50%; left:100%; transform:translateY(-50%);
    margin-left:12px; background:white; color:black; padding:10px; border-radius:6px;
    font-size:13px; line-height:1.4; width:260px; text-align:justify;
    box-shadow:0 0 10px rgba(0,0,0,0.2); z-index:10;
}
.tooltip-wrapper:hover .tooltip-bubble {display:block;}
</style>
"""

def missing_input_alert(message: str, tooltip_html: str) -> None:
    """Green boxed “please add reactants” notice with hover tooltip."""
    st.markdown(MISSING_INPUT_CSS, unsafe_allow_html=True)
    st.markdown(f"""
    <div class="tooltip-wrapper">
        <div style="
            background-color: rgba(188,212,180,0.3);
            color:#28a745; padding:0.75rem 1rem; border-radius:0.25rem;
            font-size:18px; margin-bottom:40px; box-shadow:0 2px 8px rgba(0,0,0,.15);">
            {message}
        </div>
        <div class="tooltip-bubble">{tooltip_html}</div>
    </div>
    """, unsafe_allow_html=True)

    
def dual_metric_box(
    title_left: str, txt_left: str, color_left: str, 
    title_right: str, txt_right: str, color_right: str, 
    border_color: str = "#000", title_color="#000", text_color: str = "#000") -> None:
    st.markdown(f"""
    <div style="display:flex; border:2px solid {border_color}; border-radius:10px;
                overflow:visible; margin:10px;">
      <div style="flex:1; background:{color_left}; color:{text_color}; padding:20px; text-align:center;
                  border-right:1px solid {border_color};
                  border-top-left-radius:10px; border-bottom-left-radius:10px;">
          <h4 style="margin:0; color:{title_color};">{title_left}</h4>
          <p style="margin:5px 0 0; color:{text_color};">{txt_left}</p>
      </div>
      <div style="flex:1; background:{color_right}; color:{text_color}; padding:20px; text-align:center;
                  border-top-right-radius:10px; border-bottom-right-radius:10px;">
          <h4 style="margin:0; color:{title_color};">{title_right}</h4>
          <p style="margin:5px 0 0; color:{text_color};">{txt_right}</p>
      </div>
    </div>
    """, unsafe_allow_html=True)


SOLVENT_TIP = """
• Green solvents are harmless: low-toxicity, biodegradable (e.g. water, ethanol, methanol). <br>
• Bad solvents carry harmful traits: chlorinated (toxic, ozone-depleting), <br>
  aromatic (carcinogenic), or long-chain alkanes (volatile, persistent).<br>
• Acceptable covers all others.
"""

LOGP_TIP = """
<b>Log P</b> measures molecular hydrophobicity.<br>
• Green-chemistry target: <b>1 – 3</b>.<br>
• <b>&gt; 4.5</b> ⇒ risk of persistence &amp; bio-accumulation.<br>
• <b>&lt; 1</b> ⇒ high water dispersion &amp; reduced efficacy.
"""

AA_TIP = """
• Checks reaction products for risky atoms (e.g., heavy metals, halogens). <br>
• Flags presence of concerning elements (F, Cl, Br, Pb, Hg, etc.).
"""

SAA_TIP = """
Evaluates chemical reaction products for structural concerns by checking: <br>
• Presence of long chains with more than 10 heavy atoms. <br>
• Existence or similarity to problematic chemical groups such as <br>
  carbon oxides, nitro groups, azo groups, or dichloro-aromatic groups.
"""

E_FACTOR_TIP = """
• E factor ≤ 1 → Great, almost no waste generated <br>
• E factor 1–5 → Good—manageable levels of waste typical for well-optimized processes <br>
• E factor > 5 → Poor—significant waste, process needs improvement
"""

PMI_TIP = """
• PMI < 10 → Excellent material efficiency (very little waste) <br>
• PMI 10–50 → Industry-average, acceptable but room for improvement <br>
• PMI > 50 → High waste, low material efficiency
"""

AE_M_TIP = """
• Atom Economy > 89% → ✅ Excellent incorporation of reactant mass — near-zero waste <br>
• 80–89% → ✅ Very good incorporation of reactant mass — minimal waste <br>
• 60–79% → ⚠️ Moderate incorporation of reactant mass — room for improvement <br>
• 40–59% → 🚨 Poor incorporation of reactant mass — significant waste <br>
• ≤ 39% → 🚨 Very poor incorporation of reactant mass — substantial loss <br>
"""


AE_A_TIP = """
• Atom Economy > 89% → ✅ Excellent — most atoms end up in the product <br>
• 80–89% → ✅ Very good — minimal atoms lost as by-products <br>
• 60–79% → ⚠️ Moderate — some atoms are lost as by-products <br>
• 40–59% → 🚨 Poor — many atoms wasted as by-products <br>
• ≤ 39% → 🚨 Very poor — significant atom wastage <br>
"""


def title_with_icon(title: str, tip_html: str) -> str:
    """Return `title` followed by the hoverable ℹ icon."""
    return f"{title} {tooltip_icon('ℹ', tip_html)}"