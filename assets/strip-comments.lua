-- Editing instructions in the .qmd files are written as HTML comments
-- (<!-- ... -->). This filter drops them when the site is built, so they do
-- not end up in the source of the published pages.
local function is_comment(el)
  if el.format ~= "html" then return false end
  local text = el.text:gsub("^%s+", ""):gsub("%s+$", "")
  return text:sub(1, 4) == "<!--" and text:sub(-3) == "-->" and
         not text:sub(5, -4):find("-->", 1, true)
end

function RawBlock(el)
  if is_comment(el) then return {} end
end

function RawInline(el)
  if is_comment(el) then return {} end
end
