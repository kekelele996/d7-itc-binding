# 基础镜像按补丁版本 + digest 双重锁定；换镜像时用 docker buildx imagetools inspect 取多架构 index digest。
FROM python:3.12.7-slim-bookworm@sha256:60d9996b6a8a3689d36db740b49f4327be3be09a21122bd02fb8895abb38b50d

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# 镜像构建不装任何第三方包，只用标准库；下面两处只改 apt 的取包配置：
# 软件源换成本地镜像，运行时装 agent 运行时不用再等默认源；关掉 recommends，
# 避免装 node 时被推荐包拖进整套编译工具链。判分过程不需要网络。
RUN sed -i 's|deb.debian.org|mirrors.aliyun.com|g' /etc/apt/sources.list.d/debian.sources \
    && printf 'APT::Install-Recommends "false";\nAPT::Install-Suggests "false";\n' \
        > /etc/apt/apt.conf.d/99-no-recommends

# 解题与参考解的身份由 task.toml 的 [agent] user = "runner" 指定。
# 不写 USER：Harbor 的保活进程随默认用户运行，评分前清理 runner 进程时不能把它一起结束。
RUN useradd --create-home --uid 1000 runner

WORKDIR /app
COPY app/ /app/
RUN chown -R runner:runner /app
