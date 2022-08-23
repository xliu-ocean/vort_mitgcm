%% Function to get vort
% 
function [Vor_term,Um_term_dy,Vm_term_dx] = get_vort_mitgcm_v2(Um_term,Vm_term,dxc,dyc,raz)

um_term_x = Um_term.*dxc;
vm_term_y = Vm_term.*dyc;

Um_term_dy = (um_term_x(:,2:end,:,:)-um_term_x(:,1:end-1,:,:))./raz(:,2:end);
Vm_term_dx = (vm_term_y(2:end,:,:,:)-vm_term_y(1:end-1,:,:,:))./raz(2:end,:);

% Um_term_dy = (Um_term(:,2:end,:,:)-Um_term(:,1:end-1,:,:))./dy1;
% Vm_term_dx = (Vm_term(2:end,:,:,:)-Vm_term(1:end-1,:,:,:))./dx1;
% Um_term_dy_p = 0.5*(Um_term_dy(1:end-1,:,:,:)+Um_term_dy(2:end,:,:,:));
% Vm_term_dx_p = 0.5*(Vm_term_dx(:,1:end-1,:,:)+Vm_term_dx(:,2:end,:,:));
% Um_term_dy_p = Um_term_dy(2:end,:,:,:);
% Vm_term_dx_p = Vm_term_dx(:,2:end,:,:);
Vor_term = Vm_term_dx(:,2:end,:,:)-Um_term_dy(2:end,:,:,:);
end