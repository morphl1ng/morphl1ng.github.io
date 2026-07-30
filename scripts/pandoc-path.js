/**
 * pandoc-path.js
 * 
 * Sets the correct pandoc path for Windows environments.
 * Winget installs pandoc to a user-local path that may not be in PATH.
 * This script finds pandoc automatically and configures the renderer.
 */
'use strict';
var fs = require('fs');
var path = require('path');

// Common Windows pandoc installation paths
var searchPaths = [
  'C:/Program Files/Pandoc',
  'C:/Users/' + process.env.USERNAME + '/scoop/apps/pandoc/current',
  process.env.LOCALAPPDATA + '/Microsoft/WinGet/Packages/JohnMacFarlane.Pandoc_Microsoft.Winget.Source_8wekyb3d8bbwe/pandoc-3.10',
  process.env.LOCALAPPDATA + '/Microsoft/WinGet/Packages/JohnMacFarlane.Pandoc_Microsoft.Winget.Source_8wekyb3d8bbwe/pandoc-3.9',
  process.env.LOCALAPPDATA + '/Microsoft/WinGet/Packages/JohnMacFarlane.Pandoc_Microsoft.Winget.Source_8wekyb3d8bbwe/pandoc-3.8',
];

// Also search in PATH
var pathDirs = (process.env.PATH || '').split(';');
searchPaths = searchPaths.concat(pathDirs);

var found = false;
for (var i = 0; i < searchPaths.length; i++) {
  var dir = searchPaths[i];
  var exe = path.join(dir, 'pandoc.exe');
  if (fs.existsSync(exe)) {
    hexo.config.pandoc.pandocPath = exe;
    hexo.log.info('pandoc-path: found pandoc at ' + exe);
    found = true;
    break;
  }
}

if (!found) {
  // Try just 'pandoc' - might work if in PATH
  hexo.log.warn('pandoc-path: pandoc not found via search, trying PATH');
}
