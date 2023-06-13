/*!
 * @file
 * @brief This file contains implementation of gpu
 *
 * @author Tomáš Milet, imilet@fit.vutbr.cz
 */

#include <iostream>
#include <stdlib.h>
#include <student/gpu.hpp>

typedef struct Triangle {
  OutVertex points[3];
} primitive;
struct Point {
  float x, y, z;
} point;

void clear(GPUMemory &mem, ClearCommand &cmd) {
  if (cmd.clearColor) {
    float red = cmd.color.r;
    float green = cmd.color.g;
    float blue = cmd.color.b;
    float alpha = cmd.color.a;
    for (uint32_t i = 0; i < mem.framebuffer.width * mem.framebuffer.height * mem.framebuffer.channels; i += mem.framebuffer.channels) {
      mem.framebuffer.color[i] = uint8_t(red * 255.f);
      mem.framebuffer.color[i + 1] = uint8_t(green * 255.f);
      mem.framebuffer.color[i + 2] = uint8_t(blue * 255.f);
      mem.framebuffer.color[i + 3] = uint8_t(alpha * 255.f);
    }
  }
  if (cmd.clearDepth) {
    for (uint32_t i = 0; i < mem.framebuffer.width * mem.framebuffer.height; ++i) {
      mem.framebuffer.depth[i] = cmd.depth;
    }
  }
}
void perspectiveDivision(primitive &triangle) {
  float xPos, yPos, zPos, wPos;
  for (int i = 0; i < 3; i++) {
    xPos = triangle.points[i].gl_Position.x;
    yPos = triangle.points[i].gl_Position.y;
    zPos = triangle.points[i].gl_Position.z;
    wPos = triangle.points[i].gl_Position.w;

    triangle.points[i].gl_Position.x = xPos / wPos;
    triangle.points[i].gl_Position.y = yPos / wPos;
    triangle.points[i].gl_Position.z = zPos / wPos;
  }
}
void viewportTransformation(primitive &triangle, uint32_t width, uint32_t height) {
  float xPos, yPos;
  for (int i = 0; i < 3; i++) {
    xPos = triangle.points[i].gl_Position.x;
    yPos = triangle.points[i].gl_Position.y;
    triangle.points[i].gl_Position.x = (float)width * ((xPos + 1) / 2);
    triangle.points[i].gl_Position.y = (float)height * ((yPos + 1) / 2);
  }
}
void readAttributes(GPUMemory &mem, InVertex &in, VertexArray &vao) {
  for (int i = 0; i < maxAttributes; i++) {
    const auto &attribute = vao.vertexAttrib[i];
    const void *buffer_ptr = mem.buffers[attribute.bufferID].data;
    auto buffer = buffer_ptr + attribute.offset + in.gl_VertexID * attribute.stride;
    switch (attribute.type) {
    case AttributeType::EMPTY:
      break;
    case AttributeType::FLOAT:
      in.attributes[i].v1 = *(float *)buffer;
      break;
    case AttributeType::VEC2:
      in.attributes[i].v2 = *(glm::vec2 *)buffer;
      break;
    case AttributeType::VEC3:
      in.attributes[i].v3 = *(glm::vec3 *)buffer;
      break;
    case AttributeType::VEC4:
      in.attributes[i].v4 = *(glm::vec4 *)buffer;
      break;
    case AttributeType::UINT:
      in.attributes[i].u1 = *(uint32_t *)buffer;
      break;
    case AttributeType::UVEC2:
      in.attributes[i].u2 = *(glm::uvec2 *)buffer;
      break;
    case AttributeType::UVEC3:
      in.attributes[i].u3 = *(glm::uvec3 *)buffer;
      break;
    case AttributeType::UVEC4:
      in.attributes[i].u4 = *(glm::uvec4 *)buffer;
      break;
    }
  }
}
uint32_t computeVertexID(GPUMemory &mem, VertexArray const &vao, uint32_t shaderInvocation) {
  if (vao.indexBufferID < 0)
    return shaderInvocation;

  const void *ind = mem.buffers[vao.indexBufferID].data + vao.indexOffset;

  if (vao.indexType == IndexType::UINT32) {
    return ((uint32_t *)ind)[shaderInvocation];
  } else if (vao.indexType == IndexType::UINT16) {
    return ((uint16_t *)ind)[shaderInvocation];
  } else {
    return ((uint8_t *)ind)[shaderInvocation];
  }
}
void runVertexAssembly(GPUMemory &mem, InVertex &in, VertexArray &vao, uint32_t shaderInvocation, uint32_t drawID) {
  in.gl_DrawID = drawID;
  in.gl_VertexID = computeVertexID(mem, vao, shaderInvocation);
  readAttributes(mem, in, vao);
}

void primitiveAssembly(GPUMemory &mem, uint32_t trId, VertexArray &vao, Program &prg, primitive &tr, uint32_t drawId) {
  for (int i = 0; i < 3; i++) {
    InVertex in;
    runVertexAssembly(mem, in, vao, trId * 3 + i, drawId);
    ShaderInterface si;
    si.textures = mem.textures;
    si.uniforms = mem.uniforms;
    prg.vertexShader(tr.points[i], in, si);
  }
}

float crossProduct(glm::vec2 point, glm::vec3 p1, glm::vec3 p2) { return (point.x - p2.x) * (p1.y - p2.y) - (point.y - p2.y) * (p1.x - p2.x); }

bool isPixelInTriangle(primitive triangle, float pixel_x, float pixel_y) {
  glm::vec2 point = glm::vec2(pixel_x + 0., pixel_y);
  glm::vec3 a = triangle.points[0].gl_Position;
  glm::vec3 b = triangle.points[1].gl_Position;
  glm::vec3 c = triangle.points[2].gl_Position;

  float w1 = crossProduct(point, a, b);
  float w2 = crossProduct(point, b, c);
  float w3 = crossProduct(point, c, a);

  if ((w1 > 0 && w2 > 0 && w3 > 0 || w1 < 0 && w2 < 0 && w3 < 0) && (w1 != 0 && w2 != 0 && w3 != 0)) {
    return true;
  }

  return false;
}
glm::vec4 setFrameColor(Frame &frame, primitive triangle, GPUMemory &mem, int x, int y) {
  uint8_t R, G, B, A;
  mem.framebuffer.color[(y * frame.width + x) * 4] = R;
  mem.framebuffer.color[(y * frame.width + x) * 4 + 1] = G;
  mem.framebuffer.color[(y * frame.width + x) * 4 + 2] = B;
  mem.framebuffer.color[(y * frame.width + x) * 4 + 3] = A;
  return glm::vec4(R, G, B, A);
}
int setFrameDepth(Frame &frame, int x, int y, uint8_t D) {
  frame.depth[(y * frame.width + x)] = 255;
  return frame.depth[(y * frame.width + x)] = 255;
}

void getFrameColor(Frame &frame, int *position, uint8_t color) {
  uint8_t fr_x = position[0];
  uint8_t fr_y = position[1];
  frame.color[0] = frame.color[(fr_y * frame.width + fr_x) * 4];
  frame.color[1] = frame.color[(fr_y * frame.width + fr_x) * 4 + 1];
  frame.color[2] = frame.color[(fr_y * frame.width + fr_x) * 4 + 2];
  frame.color[3] = frame.color[(fr_y * frame.width + fr_x) * 4 + 3];
}

float getFrameDepth(Frame &frame, int x, int y) {
  float depth = frame.depth[(y * frame.width + x)];
  return depth;
}

bool compare_float(float x, float y) {
  float a = 0.1f;
  if (std::abs(x - y) < a)
    return true;
  else
    return false;
}
void clampColor(OutFragment &outFragment, float min, float max) {

  if (!compare_float(outFragment.gl_FragColor.x, 1.0) && !compare_float(outFragment.gl_FragColor.y, 1.0) && !compare_float(outFragment.gl_FragColor.z, 1.0))
    outFragment.gl_FragColor.x = std::max(std::min(outFragment.gl_FragColor.x, max), min);
  outFragment.gl_FragColor.x = outFragment.gl_FragColor.x * 255;

  outFragment.gl_FragColor.y = std::max(std::min(outFragment.gl_FragColor.y, max), min);
  outFragment.gl_FragColor.y = outFragment.gl_FragColor.y * 255;

  outFragment.gl_FragColor.z = std::max(std::min(outFragment.gl_FragColor.z, max), min);
  outFragment.gl_FragColor.z = outFragment.gl_FragColor.z * 255;
}

void interpolateAttributes(primitive &primitive, Program &prg, InFragment &inFragment, float &deltaA, float &deltaB, float &deltaC) {

  // Calculate barycentric coordinates.
  float s = deltaA + deltaB + deltaC;
  float lambdaA = deltaA / s;
  float lambdaB = deltaB / s;
  float lambdaC = deltaC / s;

  // Interpolate attributes.
  for (int i = 0; i < 4; i++) {
    switch (prg.vs2fs[i]) {
    case AttributeType::EMPTY:
      break;
    case AttributeType::FLOAT:
      inFragment.attributes[i].v1 = primitive.points[0].attributes[i].v1 * lambdaA + primitive.points[1].attributes[i].v1 * lambdaB + primitive.points[2].attributes[i].v1 * lambdaC;
      break;
    case AttributeType::VEC2:
      inFragment.attributes[i].v2 = primitive.points[0].attributes[i].v2 * lambdaA + primitive.points[1].attributes[i].v2 * lambdaB + primitive.points[2].attributes[i].v2 * lambdaC;
      break;
    case AttributeType::VEC3:
      inFragment.attributes[i].v3 = primitive.points[0].attributes[i].v3 * lambdaA + primitive.points[1].attributes[i].v3 * lambdaB + primitive.points[2].attributes[i].v3 * lambdaC;
      break;
    case AttributeType::VEC4:
      inFragment.attributes[i].v4 = primitive.points[0].attributes[i].v4 * lambdaA + primitive.points[1].attributes[i].v4 * lambdaB + primitive.points[2].attributes[i].v4 * lambdaC;
      break;
    }
  }

  // Calculate fragment depth.
  float z0 = primitive.points[0].gl_Position.z;
  float z1 = primitive.points[1].gl_Position.z;
  float z2 = primitive.points[2].gl_Position.z;
  inFragment.gl_FragCoord.z = z0 * lambdaA + z1 * lambdaB + z2 * lambdaC;
}

glm::vec3 inTriangleBarycentric(primitive &primitive, float x, float y) {
  glm::vec4 p0 = primitive.points[0].gl_Position;
  glm::vec4 p1 = primitive.points[1].gl_Position;
  glm::vec4 p2 = primitive.points[2].gl_Position;
  float lambda0 = ((p1.y - p2.y) * (x - p2.x) + (p2.x - p1.x) * (y - p2.y)) / ((p1.y - p2.y) * (p0.x - p2.x) + (p2.x - p1.x) * (p0.y - p2.y));
  float lambda1 = ((p2.y - p0.y) * (x - p2.x) + (p0.x - p2.x) * (y - p2.y)) / ((p1.y - p2.y) * (p0.x - p2.x) + (p2.x - p1.x) * (p0.y - p2.y));
  float lambda2 = 1.0f - lambda0 - lambda1;
  float z = p0.z * lambda0 + p1.z * lambda1 + p2.z * lambda2;
  return glm::vec3(lambda0, lambda1, lambda2);
}

void getEdges(primitive &primitive, glm::vec3 edges[]) {
  edges[0] = {primitive.points[1].gl_Position.x - primitive.points[0].gl_Position.x, primitive.points[1].gl_Position.y - primitive.points[0].gl_Position.y, 0};
  edges[1] = {primitive.points[2].gl_Position.x - primitive.points[1].gl_Position.x, primitive.points[2].gl_Position.y - primitive.points[1].gl_Position.y, 0};
  edges[2] = {primitive.points[0].gl_Position.x - primitive.points[2].gl_Position.x, primitive.points[0].gl_Position.y - primitive.points[2].gl_Position.y, 0};
}
void createFragment(GPUMemory &mem, primitive &primitive, InFragment &in, float x, float y) {
  glm::vec3 edges[3];
  getEdges(primitive, edges);
  in.gl_FragCoord = glm::vec4(x, y, 0, 0);
  auto deltas = inTriangleBarycentric(primitive, x, y);
  in.gl_FragCoord.z = deltas[0] * primitive.points[0].gl_Position.z + deltas[1] * primitive.points[1].gl_Position.z + deltas[2] * primitive.points[2].gl_Position.z;
  return;
}

void setFragmentDepth(glm::vec4 &fragment_color, Triangle triangle, float lambda, float lambda_x, float lambda_y, float lambda_z, Frame &frame, int x, int y) {
  glm::vec3 A_color = triangle.points[0].attributes[0].v3;
  glm::vec3 B_color = triangle.points[1].attributes[0].v3;
  glm::vec3 C_color = triangle.points[2].attributes[0].v3;
  float h0 = triangle.points[0].gl_Position.w;
  float h1 = triangle.points[1].gl_Position.w;
  float h2 = triangle.points[2].gl_Position.w;
  fragment_color = glm::vec4((A_color * lambda_x / h0 + B_color * lambda_y / h1 + C_color * lambda_z / h2) / (lambda / h0 + lambda / h1 + lambda / h2), 1);
}

bool counterClockWise(glm::vec3 edges[]) {
  if (glm::cross(edges[0], edges[1]).z >= 0 && glm::cross(edges[1], edges[2]).z >= 0 && glm::cross(edges[2], edges[0]).z >= 0)
    return true;
  else
    return false;
}
void perFragmentOperations(Frame &frame, uint32_t x, uint32_t y, OutFragment &outFragment, float depth) {
  uint32_t pos = (frame.height) * y + x;
  if (frame.depth[pos] > depth) {
    if (outFragment.gl_FragColor.w > 0.5f) {
      frame.depth[pos] = depth;
    }
    frame.color[4 * pos] = frame.color[4 * pos] * (1.0f - outFragment.gl_FragColor.w) + outFragment.gl_FragColor.x * outFragment.gl_FragColor.w;
    frame.color[4 * pos + 1] = frame.color[4 * pos + 1] * (1.0f - outFragment.gl_FragColor.w) + outFragment.gl_FragColor.y * outFragment.gl_FragColor.w;
    frame.color[4 * pos + 2] = frame.color[4 * pos + 2] * (1.0f - outFragment.gl_FragColor.w) + outFragment.gl_FragColor.z * outFragment.gl_FragColor.w;
  }
}
void rasterizeTriangle(DrawCommand &cmd, GPUMemory &mem, primitive &triangle, Program &prg) {
  glm::vec3 edges[3];
  getEdges(triangle, edges);
  for (int x = 0; x < mem.framebuffer.width; x++) {
    for (int y = 0; y < mem.framebuffer.height; y++) {
      float pointX = (float)x + 0.5;
      float pointY = (float)y + 0.5;
      if (isPixelInTriangle(triangle, pointX, pointY)) {
        if (!cmd.backfaceCulling || counterClockWise(edges)) {
          InFragment inFragment;
          OutFragment outFragment;
          Frame frame = mem.framebuffer;
          createFragment(mem, triangle, inFragment, pointX, pointY);
          auto deltas = inTriangleBarycentric(triangle, x, y);
          interpolateAttributes(triangle, prg, inFragment, deltas[0], deltas[1], deltas[2]);
          float depth = getFrameDepth(mem.framebuffer, x, y);
          float fragDepth = setFrameDepth(frame, x, y, inFragment.gl_FragCoord.z);
          auto fragColor = setFrameColor(frame, triangle, mem, x, y);
          clampColor(outFragment, 0, 1);
          perFragmentOperations(frame, x, y, outFragment, inFragment.gl_FragCoord.z);
          ShaderInterface si;
          si.textures = mem.textures;
          si.uniforms = mem.uniforms;

          prg.fragmentShader(outFragment, inFragment, si);
        }
      }
    }
  }
}

void draw(GPUMemory &mem, DrawCommand &cmd, uint32_t drawID) {
  Program prg = mem.programs[cmd.programID];
  VertexShader vs = prg.vertexShader;
  VertexArray vao = cmd.vao;

  // Iterate over all vertices
  for (int trId = 0; trId < cmd.nofVertices / 3; trId++) {
    primitive tr;
    primitiveAssembly(mem, trId, vao, prg, tr, drawID);
    perspectiveDivision(tr);
    viewportTransformation(tr, mem.framebuffer.width, mem.framebuffer.height);
    rasterizeTriangle(cmd, mem, tr, prg);
  }
}

void loadVertex(GPUMemory &mem, InVertex &inVertex, VertexArray vao, uint32_t ID, bool userData, uint8_t *n_data) {
  uint64_t offset;
  uint64_t stride;
  AttributeType type;
  uint8_t *data;
  for (int j = 0; j < maxAttributes; j++) {
    type = vao.vertexAttrib[j].type;
    offset = vao.vertexAttrib[j].offset;
    stride = vao.vertexAttrib[j].stride;
    if (userData)
      data = (uint8_t *)vao.vertexAttrib[j].bufferID;
    else
      data = (uint8_t *)n_data;
    switch (type) {
    case AttributeType::FLOAT:
      inVertex.attributes[j].v1 = (float)((float *)(data + offset + inVertex.gl_VertexID * stride))[0];
      break;
    case AttributeType::VEC2:
      inVertex.attributes[j].v2 = (glm::vec2)((glm::vec2 *)(data + offset + inVertex.gl_VertexID * stride))[0];
      break;
    case AttributeType::VEC3:
      inVertex.attributes[j].v3 = (glm::vec3)((glm::vec3 *)(data + offset + inVertex.gl_VertexID * stride))[0];
      break;
    case AttributeType::VEC4:
      inVertex.attributes[j].v4 = (glm::vec4)((glm::vec4 *)(data + offset + inVertex.gl_VertexID * stride))[0];
      break;
    default:
      continue;
    }
  }
  return;
}

//! [gpu_execute]
void gpu_execute(GPUMemory &mem, CommandBuffer &cb) {
  (void)mem;
  (void)cb;
  int draw_num = 0;
  for (uint32_t i = 0; i < cb.nofCommands; ++i) {
    CommandType type = cb.commands[i].type;
    auto data = cb.commands[i].data;

    if (type == CommandType::DRAW) {
      draw(mem, data.drawCommand, draw_num);
      draw_num++;
    }
    if (type == CommandType::CLEAR) {
      clear(mem, data.clearCommand);
    }
  }
}
/// \todo Tato funkce reprezentuje funkcionalitu grafické karty.<br>
/// Měla by umět zpracovat command buffer, čistit framebuffer a kresli.<br>
/// mem obsahuje paměť grafické karty.
/// cb obsahuje command buffer pro zpracování.
/// Bližší informace jsou uvedeny na hlavní stránce dokumentace.

//! [gpu_execute]

/**
 * @brief This function reads color from texture.
 *
 * @param texture texture
 * @param uv uv coordinates
 *
 * @return color 4 floats
 */
glm::vec4 read_texture(Texture const &texture, glm::vec2 uv) {
  if (!texture.data)
    return glm::vec4(0.f);
  auto uv1 = glm::fract(uv);
  auto uv2 = uv1 * glm::vec2(texture.width - 1, texture.height - 1) + 0.5f;
  auto pix = glm::uvec2(uv2);
  // auto t   = glm::fract(uv2);
  glm::vec4 color = glm::vec4(0.f, 0.f, 0.f, 1.f);
  for (uint32_t c = 0; c < texture.channels; ++c)
    color[c] = texture.data[(pix.y * texture.width + pix.x) * texture.channels + c] / 255.f;
  return color;
}
