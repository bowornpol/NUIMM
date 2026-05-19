files <- list.files("d:/NUIMM/NUIMM_web", pattern = "\\.html$", full.names = TRUE)

for (f in files) {
  content <- readLines(f, warn = FALSE)
  content_str <- paste(content, collapse = "\n")
  
  # Remove the footer-right block completely (supports both single line and multi-line formats)
  content_str <- gsub('<div class="footer-right">.*?</div></div></div></div>', '', content_str)
  content_str <- gsub('<div class="footer-right">\\s*<div class="brand-pill">.*?</div\\s*>\\s*</div\\s*>\\s*</div\\s*>', '', content_str)
  # Actually, let's just do a greedy replace from <div class="footer-right"> up to Mahidol University</div></div></div></div>
  content_str <- sub('<div class="footer-right">.*Mahidol University</div>\\s*</div>\\s*</div>\\s*</div>', '', content_str)
  
  # Update Header Logo for non-viewer files
  if (!grepl("NUIMM_viewer.html", f)) {
    content_str <- sub(
      '<a href="NUIMM_homepage.html" class="logo-container">\\s*<div class="graphic-b">.*?<div class="logo-text">NUIMM</div>\\s*</a>',
      '<div style="display: flex; align-items: center; gap: 20px;">\n                <a href="NUIMM_homepage.html" class="logo-container">\n                    <div class="graphic-b">\n                        <div class="net-line nl1"></div>\n                        <div class="net-line nl2"></div>\n                        <div class="net-line nl3"></div>\n                        <div class="net-line nl4"></div>\n                        <div class="sphere-node sn1"></div>\n                        <div class="sphere-node sn2"></div>\n                        <div class="sphere-node sn3"></div>\n                    </div>\n                    <div class="logo-text">NUIMM</div>\n                </a>\n                <img src="logo_mu.png" alt="Mahidol University" style="height: 38px;">\n            </div>',
      content_str
    )
  } else {
    # Update Header Logo for viewer file
    content_str <- sub(
      '<a href="NUIMM_homepage.html" class="logo-container">\\s*<div class="graphic-b">.*?<div class="logo-text">NUIMM Result Viewer</div>\\s*</a>',
      '<div style="display: flex; align-items: center; gap: 20px;">\n            <a href="NUIMM_homepage.html" class="logo-container">\n                <div class="graphic-b">\n                    <div class="net-line nl1"></div>\n                    <div class="net-line nl2"></div>\n                    <div class="net-line nl3"></div>\n                    <div class="sphere-node sn1"></div>\n                    <div class="sphere-node sn2"></div>\n                    <div class="sphere-node sn3"></div>\n                </div>\n                <div class="logo-text">NUIMM Result Viewer</div>\n            </a>\n            <img src="logo_mu.png" alt="Mahidol University" style="height: 38px;">\n        </div>',
      content_str
    )
  }
  
  writeLines(content_str, f)
}

print("Update complete")
