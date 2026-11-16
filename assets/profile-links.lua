-- arXiv and ORCID links, shown only once their ids are filled in.
-- The ids are the two "profile-links" settings in _quarto.yml. While a setting is
-- empty its link is not written anywhere. Once it is filled in, the link appears
--   * after the Email and GitHub buttons, in every block marked {.wm-profile-links}
--     (the home page and the CV page), and
--   * next to the Email and GitHub icons in the menu bar.
local NAME = "William Mahaney"
local SITES = {
  { key = "arxiv-author", text = "arXiv", url = "https://arxiv.org/a/", label = NAME .. " on arXiv" },
  { key = "orcid", text = "ORCID", url = "https://orcid.org/", label = NAME .. " on ORCID" },
}
local links = {}

local function read_settings(meta)
  local cfg = meta["profile-links"]
  if type(cfg) ~= "table" then return end
  for _, site in ipairs(SITES) do
    local id = cfg[site.key] and pandoc.utils.stringify(cfg[site.key]) or ""
    id = id:gsub("^%s+", ""):gsub("%s+$", "")
    if id ~= "" then
      if id:match("^[%w_%-]+$") then
        links[#links + 1] = { href = site.url .. id, text = site.text, label = site.label }
      else
        io.stderr:write("profile-links: ignoring " .. site.key .. " = '" .. id ..
          "' (letters, digits, '-' and '_' only; give the id, not the full address)\n")
      end
    end
  end
end

-- The menu bar is not part of the page text, so its links are added by a few lines of
-- script, written into the page only when there is a link to add.
local function navbar_script()
  local items = {}
  for _, l in ipairs(links) do
    items[#items + 1] = string.format('["%s","%s","%s"]', l.href, l.label, l.text)
  end
  return '<script>document.addEventListener("DOMContentLoaded",function(){' ..
    'var bar=document.querySelector(".quarto-navbar-tools");if(!bar)return;' ..
    'var before=bar.querySelector(".quarto-color-scheme-toggle");' ..
    '[' .. table.concat(items, ",") .. '].forEach(function(l){' ..
    'var a=document.createElement("a");a.href=l[0];a.className="quarto-navigation-tool wm-nav-id px-1";' ..
    'a.setAttribute("aria-label",l[1]);a.title=l[1];a.textContent=l[2];bar.insertBefore(a,before);});});</script>'
end

return {
  {
    Meta = function(meta)
      links = {}
      read_settings(meta)
      if #links == 0 then return nil end
      local includes = meta["header-includes"]
      local list = pandoc.List({})
      if includes ~= nil then
        if pandoc.utils.type(includes) == "List" then list:extend(includes) else list:insert(includes) end
      end
      list:insert(pandoc.RawBlock("html", navbar_script()))
      meta["header-includes"] = list
      return meta
    end,
  },
  {
    Div = function(div)
      if #links == 0 or not div.classes:includes("wm-profile-links") then return nil end
      local buttons = {}
      for _, l in ipairs(links) do
        buttons[#buttons + 1] = pandoc.Space()
        buttons[#buttons + 1] = pandoc.Link(l.text, l.href, l.label,
          pandoc.Attr("", { "wm-btn", "wm-btn-outline" }, { ["aria-label"] = l.label }))
      end
      local last = div.content[#div.content]
      if last and (last.t == "Para" or last.t == "Plain") then
        last.content:extend(buttons)
      else
        table.remove(buttons, 1)
        div.content:insert(pandoc.Para(buttons))
      end
      return div
    end,
  },
}
