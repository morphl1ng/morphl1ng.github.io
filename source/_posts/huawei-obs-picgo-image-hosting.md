---
title: 用华为云 OBS + PicGo 搭建个人图床
date: 2026-07-23
tags:
  - 图床
  - PicGo
  - 华为云
categories:
  - 技术笔记
---

记录图床构建

## 一、华为云 OBS 配置

### 1. 创建桶

1. 登录 [对象存储服务OBS_官网_云存储服务_数据云存储解决方案-华为云](https://www.huaweicloud.com/product/obs.html)，完成实名认证后进入 **对象存储服务 OBS**

   ![image-20260723142015405](https://morphl1ng-blog.obs.cn-east-3.myhuaweicloud.com/blog/image-20260723142015405.png)

   ![image-20260723142104566](https://morphl1ng-blog.obs.cn-east-3.myhuaweicloud.com/blog/image-20260723142104566.png)

2. 点击 **创建桶**，关键选项：
   - **区域**：就近选，例如 `华东-上海一`

   - **存储类别**：标准存储

   - **桶名称**：全局唯一，例如 `myblog-你的id`

   - **开启公共读**

     ![image-20260723142444349](https://morphl1ng-blog.obs.cn-east-3.myhuaweicloud.com/blog/image-20260723142444349.png)

### 2. 开启公共读（关键！）

图床图片必须让所有人能通过链接访问，否则外链全部 403：

进入桶 → **权限控制** → **桶策略** → 创建桶策略设置为 **公共读**（匿名用户可读取桶内对象）。

![image-20260723142706307](https://morphl1ng-blog.obs.cn-east-3.myhuaweicloud.com/blog/image-20260723142706307.png)

### 3. 获取访问密钥 AK/SK、查看Endpoint

右上角账号名 → **我的凭证** → **访问密钥** → **新增访问密钥**，保存好下载的 `credentials.csv`。

> ⚠️ 安全建议：AK/SK 等于账号密码，**永远不要提交进 git 仓库**。更稳妥的做法是在 IAM 里创建一个子用户，只授予 OBS 的读写权限，用子用户的密钥。

在概览中查看Endpoint，示例（华东-上海一）obs.cn-east-3.myhuaweicloud.com

## 二、PicGo 配置

### 1. 安装

从 [PicGo Releases](https://github.com/Molunerfinn/PicGo/releases) 下载 Windows 安装包安装。

### 2. 安装华为云上传插件

PicGo 本体不带华为云支持，进入 **插件设置**，搜索 `huawei`，安装华为云 OBS 专用插件。

![image-20260723143037314](https://morphl1ng-blog.obs.cn-east-3.myhuaweicloud.com/blog/image-20260723143037314.png)

### 3. 配置华为云OBS

![image-20260723143148886](https://morphl1ng-blog.obs.cn-east-3.myhuaweicloud.com/blog/image-20260723143148886.png)

## 三、Typora 配置

修改好后右键点击上传图像，链接会自动切换。

![image-20260723143358067](https://morphl1ng-blog.obs.cn-east-3.myhuaweicloud.com/blog/image-20260723143358067.png)

