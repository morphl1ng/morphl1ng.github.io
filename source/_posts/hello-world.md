---
title: 博客开张：它是怎么搭起来的
date: 2026-07-23
tags:
  - 博客
  - Hexo
categories:
  - 随笔
---

这是我的第一篇文章。这个博客基于 **Hexo + NexT 主题 + GitHub Pages**，由 GitHub Actions 自动部署。

## 日常写作流程

```bash
# 新建一篇文章
hexo new post "文章标题"

# 本地预览（http://localhost:4000）
hexo server

# 写完直接推送，GitHub Actions 会自动构建并上线
git add .
git commit -m "新文章"
git push
```

几分钟后刷新 https://morphl1ng.github.io 就能看到新文章。

## 文章头部模板

每篇文章开头的 Front-matter 控制标题、日期、分类标签：

```yaml
---
title: 文章标题
date: 2026-07-23
tags:
  - 标签1
  - 标签2
categories:
  - 分类名
---
```

## 常用目录

| 位置 | 作用 |
| --- | --- |
| `source/_posts/` | 所有文章 |
| `_config.yml` | 站点主配置 |
| `_config.next.yml` | 主题配置 |
| `source/about/index.md` | 关于页 |

开始写作吧！
