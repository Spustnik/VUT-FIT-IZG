/*!
 * @file
 * @brief This file contains functions for model rendering
 *
 * @author Tomáš Milet, imilet@fit.vutbr.cz
 */
#include <student/drawModel.hpp>
#include <student/gpu.hpp>

///\endcond
void addDrawCommand(CommandBuffer &cb, size_t nofVertices, Mesh mesh) {
  (void)cb;
  cb.commands[cb.nofCommands].type = CommandType::DRAW;
  cb.commands[cb.nofCommands].data.drawCommand.nofVertices = nofVertices;
  cb.commands[cb.nofCommands].data.drawCommand.backfaceCulling = true;
  cb.commands[cb.nofCommands].data.drawCommand.programID = 0;
  cb.nofCommands++;
  if (mesh.doubleSided) {
    cb.commands[cb.nofCommands].data.drawCommand.backfaceCulling = false;
  }
}

void prepareNode(GPUMemory &mem, CommandBuffer &cb, Model const &model, Node const &node, glm::mat4 const &parentMatrix) {
  if (node.mesh >= 0) {
    Mesh const &mesh = model.meshes[node.mesh];
    glm::mat4 modelMatrix = parentMatrix * node.modelMatrix;
    mem.uniforms[10 + ((cb.nofCommands - 1) * 5)].m4 = modelMatrix;
    modelMatrix = glm::transpose(glm::inverse(parentMatrix * node.modelMatrix));
    mem.uniforms[10 + ((cb.nofCommands - 1) * 5) + 1].m4 = modelMatrix;
    glm::vec4 modelMatrix2 = glm::vec4(model.meshes[node.mesh].diffuseColor.r, model.meshes[node.mesh].diffuseColor.g, model.meshes[node.mesh].diffuseColor.b, model.meshes[node.mesh].diffuseColor.a);
    mem.uniforms[10 + ((cb.nofCommands - 1) * 5) + 2].v4 = modelMatrix2;
    auto modelMatrix3 = glm::int32_t(model.meshes[node.mesh].diffuseTexture);
    mem.uniforms[10 + ((cb.nofCommands - 1) * 5) + 3].i1 = modelMatrix3;
    auto modelMatrix4 = glm::float32(model.meshes[node.mesh].doubleSided);
    mem.uniforms[10 + ((cb.nofCommands - 1) * 5) + 4].v1 = modelMatrix4;
    addDrawCommand(cb, mesh.nofIndices, mesh);
  }
  for (size_t i = 0; i < node.children.size(); ++i) {
    prepareNode(mem, cb, model, node.children[i], parentMatrix * node.modelMatrix);
  }
}
/**
 * @brief This function prepares model into memory and creates command buffer
 *
 * @param mem gpu memory
 * @param commandBuffer command buffer
 * @param model model structure
 */
//! [drawModel]
void prepareModel(GPUMemory &mem, CommandBuffer &cb, Model const &model) {
  (void)mem;
  (void)cb;
  (void)model;
  for (int i = 0; i < model.buffers.size(); i++) {
    mem.buffers[i].data = model.buffers[i].data;
    mem.buffers[i].size = model.buffers[i].size;
  }
  for (int i = 0; i < model.textures.size(); i++) {
    mem.textures[i].data = model.textures[i].data;
    mem.textures[i].width = model.textures[i].width;
    mem.textures[i].height = model.textures[i].height;
    mem.textures[i].channels = model.textures[i].channels;
  }
  mem.programs[0].vs2fs[0] = AttributeType::VEC3;
  mem.programs[0].vs2fs[1] = AttributeType::VEC3;
  mem.programs[0].vs2fs[2] = AttributeType::VEC2;
  mem.programs[0].vs2fs[3] = AttributeType::UINT;
  mem.programs[0].vertexShader = drawModel_vertexShader;
  mem.programs[0].fragmentShader = drawModel_fragmentShader;
  cb.commands[0].type = CommandType::CLEAR;
  cb.commands[0].data.clearCommand.clearColor = true;
  cb.commands[0].data.clearCommand.clearDepth = true;
  cb.commands[0].data.clearCommand.color = glm::vec4(0.1, 0.15, 0.1, 1);
  cb.commands[0].data.clearCommand.depth = 1e+11f;
  cb.nofCommands = 1;

  for (size_t i = 0; i < model.roots.size(); ++i) {
    prepareNode(mem, cb, model, model.roots[i], glm::mat4(1.f));
  }
}
/// \todo Tato funkce připraví command buffer pro model a nastaví správně pamět grafické karty.<br>
/// Vaším úkolem je správně projít model a vložit vykreslovací příkazy do commandBufferu.
/// Zároveň musíte vložit do paměti textury, buffery a uniformní proměnné, které buffer command buffer využívat.
/// Bližší informace jsou uvedeny na hlavní stránce dokumentace a v testech.
//! [drawModel]

/**
 * @brief This function represents vertex shader of texture rendering method.
 *
 * @param outVertex output vertex
 * @param inVertex input vertex
 * @param si shader interface
 */
//! [drawModel_vs]
void drawModel_vertexShader(OutVertex &outVertex, InVertex const &inVertex, ShaderInterface const &si) {
  (void)outVertex;
  (void)inVertex;
  (void)si;

  /// \todo Tato funkce reprezentujte vertex shader.<br>
  /// Vaším úkolem je správně trasnformovat vrcholy modelu.
  /// Bližší informace jsou uvedeny na hlavní stránce dokumentace.
}
//! [drawModel_vs]

/**
 * @brief This functionrepresents fragment shader of texture rendering method.
 *
 * @param outFragment output fragment
 * @param inFragment input fragment
 * @param si shader interface
 */
//! [drawModel_fs]
void drawModel_fragmentShader(OutFragment &outFragment, InFragment const &inFragment, ShaderInterface const &si) {
  (void)outFragment;
  (void)inFragment;
  (void)si;

  // outFragment.gl_FragColor = glm::vec4(glm::vec3(aL + dL), diffuseColor.a);
  /// \todo Tato funkce reprezentujte fragment shader.<br>
  /// Vaším úkolem je správně obarvit fragmenty a osvětlit je pomocí lambertova osvětlovacího modelu.
  /// Bližší informace jsou uvedeny na hlavní stránce dokumentace.
}

//! [drawModel_fs]
